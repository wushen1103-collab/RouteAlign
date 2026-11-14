#!/usr/bin/env python3
"""Endpoint-only evaluation of PB-valid receptor/consensus group selections."""
from __future__ import annotations
from routealign import settings

import argparse
import atexit
import csv
import hashlib
import json
import os
import pathlib
import shutil
import signal
import tempfile
import time

from rdkit import Chem, rdBase
from posebusters import __version__ as posebusters_version
from posebusters.modules.rmsd import check_rmsd

ROOT = settings.workspace()
MANIFEST = ROOT / "config/cohort.csv"
DATASET = ROOT / "data/AlignDockBench"
SELECTOR = ROOT / "runs/receptor_select"
REPAIR = ROOT / "runs/receptor_repair"
PREREG = ROOT / "config/protocol.json"
OUTPUT = ROOT / "runs/receptor_evaluate"
EXPECTED_MANIFEST_SHA256 = "4d9924a2985c556ca45060a7201d6f21bce91021d33c53206dd8f048662c9549"
ACTIVE: pathlib.Path | None = None

def sha256(path: pathlib.Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def atomic_json(path: pathlib.Path, payload: dict) -> None:
    with path.open('w',encoding='utf-8') as f:
        json.dump(payload,f,indent=2,sort_keys=True,allow_nan=False); f.write('\n'); f.flush(); os.fsync(f.fileno())

def clean() -> None:
    global ACTIVE
    if ACTIVE is not None and ACTIVE.exists(): shutil.rmtree(ACTIVE)
    ACTIVE=None

def stop(sig,_frame) -> None: clean(); raise SystemExit(128+sig)
atexit.register(clean); signal.signal(signal.SIGTERM,stop); signal.signal(signal.SIGINT,stop)

def task_table():
    if sha256(MANIFEST)!=EXPECTED_MANIFEST_SHA256: raise RuntimeError('manifest identity drift')
    prereg=json.loads(PREREG.read_text(encoding='utf-8'))
    pairs=[x for x in prereg['pairs'] if not x['hard_low_2d']]; groups=[tuple(x) for x in prereg['outer_seed_groups']]
    with MANIFEST.open(newline='',encoding='utf-8') as f: rows=[x for x in csv.DictReader(f) if x['final_primary_eligible']=='True']
    by={f"{x['family']}__{x['query_id']}":x for x in rows}
    if len(pairs)!=18 or len(groups)!=5: raise RuntimeError('frozen non-low2d design mismatch')
    return [(i,p,gi,g,by[p['pair']]) for i,(p,gi,g) in enumerate((p,gi,g) for p in pairs for gi,g in enumerate(groups,1))]

def load_single(path):
    ms=[m for m in Chem.SDMolSupplier(str(path),removeHs=True,sanitize=True) if m is not None]
    if len(ms)!=1 or ms[0].GetNumConformers()!=1: raise RuntimeError('pose parse failed: '+str(path))
    return ms[0]

def run(index,audit_only=False,fault_after_metric=False,hold_after_metric=0.0):
    global ACTIVE
    table=task_table()
    if pathlib.Path.cwd().resolve()!=ROOT or not 0<=index<len(table): raise RuntimeError('scope/index mismatch')
    _,pair_spec,group_index,seeds,row=table[index]; pair=pair_spec['pair']
    pose=SELECTOR/pair/f'group_{group_index}'/'selected.sdf'; selection_path=SELECTOR/pair/f'group_{group_index}'/'result.json'; native=DATASET/row['family']/row['query_id']/'ligand.sdf'
    for p in (pose,selection_path,native):
        p.resolve().relative_to(ROOT)
        if not p.is_file() or p.is_symlink(): raise FileNotFoundError(p)
    selection=json.loads(selection_path.read_text(encoding='utf-8')); selected_seed=int(selection['selected_seed'])
    if selected_seed not in seeds: raise RuntimeError('selected seed outside frozen group')
    repair_path=REPAIR/pair/f'seed_{selected_seed}'/'metadata.json'; repair=json.loads(repair_path.read_text(encoding='utf-8'))
    expected={'task_index':index,'pair':pair,'group_index':group_index,'selected_seed':selected_seed,'hard_low_2d':False,'hard_flex':bool(pair_spec['hard_flex']),'hard_scaffold':bool(pair_spec['hard_scaffold']),'manifest_sha256':EXPECTED_MANIFEST_SHA256,'preregister_sha256':sha256(PREREG),'prediction_sha256':sha256(pose),'selector_result_sha256':sha256(selection_path),'repair_metadata_sha256':sha256(repair_path),'native_sha256':sha256(native),'worker_sha256':sha256(pathlib.Path(__file__))}
    if audit_only:return {**expected,'status':'audit_pass','endpoint_only':True,'threads_max':1,'gpu_hidden':True}
    if os.environ.get('CUDA_VISIBLE_DEVICES')!='' or os.getpid()!=os.getpgrp(): raise RuntimeError('GPU or process-group isolation missing')
    task_root=OUTPUT/pair/f'group_{group_index}'; result_path=task_root/'result.json'
    if result_path.is_file():
        payload=json.loads(result_path.read_text(encoding='utf-8'))
        if all(payload.get(k)==v for k,v in expected.items()): return {**expected,'status':'resumed','rmsd_A':payload['rmsd_A']}
        raise RuntimeError('resume identity drift')
    metric=check_rmsd(load_single(pose),load_single(native),rmsd_threshold=2.0,heavy_only=True,choose_by='rmsd')['results']; rmsd=float(metric['rmsd'])
    if fault_after_metric: raise SystemExit(97)
    if hold_after_metric>0: time.sleep(hold_after_metric)
    task_root.parent.mkdir(parents=True,exist_ok=True); stage=pathlib.Path(tempfile.mkdtemp(prefix=task_root.name+'.staging.',dir=task_root.parent)); ACTIVE=stage
    success=bool(metric['rmsd_within_threshold']); pb=bool(repair['final_PB_pass'])
    payload={**expected,'status':'success','rmsd_A':rmsd,'success_le_2A':success,'PB20_pass':pb,'joint_success':bool(success and pb),'repair_decision':repair['decision'],'query_native_ligand_read_for_endpoint':True,'selector_read_native':False,'posebusters_version':posebusters_version,'rdkit_version':rdBase.rdkitVersion,'SOTA':False}
    atomic_json(stage/'result.json',payload); os.replace(stage,task_root); ACTIVE=None
    return {**expected,'status':'computed','rmsd_A':rmsd,'success_le_2A':success,'PB20_pass':pb,'joint_success':bool(success and pb)}

def main():
    p=argparse.ArgumentParser(); p.add_argument('--index',type=int,required=True); p.add_argument('--audit-only',action='store_true'); p.add_argument('--fault-after-metric',action='store_true'); p.add_argument('--hold-after-metric',type=float,default=0.0); a=p.parse_args(); print(json.dumps(run(a.index,a.audit_only,a.fault_after_metric,a.hold_after_metric),sort_keys=True))
if __name__=='__main__':main()
