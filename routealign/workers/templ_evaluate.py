#!/usr/bin/env python3
"""Symmetry-aware RMSD evaluation of the same-row TEMPL author-core rerun."""
from __future__ import annotations
from routealign import settings
import argparse,csv,hashlib,json,os,pathlib,shutil,signal,tempfile,time,atexit
from rdkit import Chem,rdBase
from posebusters import __version__ as posebusters_version
from posebusters.modules.rmsd import check_rmsd

ROOT = settings.workspace()
MANIFEST=ROOT/'config/cohort.csv'
PRED=ROOT/'runs/templ_generate'
DATA=ROOT/'data/AlignDockBench'
OUT=pathlib.Path(os.environ.get('ROUTEALIGN_ADB_TEMPL_EVAL_OUTPUT',str(ROOT/'runs/templ_evaluate')))
MSH='4d9924a2985c556ca45060a7201d6f21bce91021d33c53206dd8f048662c9549'
ACTIVE=None
def sha(p):
 h=hashlib.sha256()
 with pathlib.Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def atomic(p,d):
 p.parent.mkdir(parents=True,exist_ok=True);fd,t=tempfile.mkstemp(prefix=p.name+'.',suffix='.tmp',dir=p.parent)
 try:
  with os.fdopen(fd,'w') as f:json.dump(d,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
  os.replace(t,p)
 finally:
  if os.path.exists(t):os.unlink(t)
def clean():
 global ACTIVE
 if ACTIVE and ACTIVE.exists():shutil.rmtree(ACTIVE,ignore_errors=True)
 ACTIVE=None
def stop(sig,_):clean();raise SystemExit(128+sig)
atexit.register(clean);signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
def rows():
 if sha(MANIFEST)!=MSH:raise RuntimeError('manifest drift')
 with MANIFEST.open(newline='') as f:r=[x for x in csv.DictReader(f) if x['final_primary_eligible']=='True']
 r.sort(key=lambda x:(x['family'],x['query_id']))
 if len(r)!=44:raise RuntimeError('cohort')
 return r
def mol(path):
 values=[m for m in Chem.SDMolSupplier(str(path),removeHs=True,sanitize=True) if m is not None]
 if len(values)!=1 or values[0].GetNumConformers()!=1:raise RuntimeError('one conformer required '+str(path))
 return values[0]
def run(i,fault=False,hold=0.):
 global ACTIVE
 table=rows()
 if pathlib.Path.cwd().resolve()!=ROOT or not 0<=i<len(table) or (os.cpu_count() or 0)<31:raise RuntimeError('scope/resources')
 r=table[i];name=f"{r['family']}__{r['query_id']}";pred=PRED/name/'selected.sdf';native=DATA/r['family']/r['query_id']/'ligand.sdf';meta=PRED/name/'result.json'
 for p in (pred,native,meta):
  p.resolve().relative_to(ROOT); 
  if not p.is_file():raise FileNotFoundError(p)
 expected={'task_index':i,'family':r['family'],'query_id':r['query_id'],'manifest_sha256':MSH,'prediction_sha256':sha(pred),'native_sha256':sha(native),'generation_result_sha256':sha(meta),'worker_sha256':sha(__file__)}
 task=OUT/name;rp=task/'result.json'
 if rp.is_file():
  d=json.loads(rp.read_text())
  if all(d.get(k)==v for k,v in expected.items()):return {**expected,'status':'resumed','rmsd_A':d['rmsd_A']}
  raise RuntimeError('resume drift')
 checked=check_rmsd(mol(pred),mol(native),rmsd_threshold=2.0,heavy_only=True,choose_by='rmsd')['results']
 if fault:raise SystemExit(97)
 if hold:time.sleep(hold)
 task.parent.mkdir(parents=True,exist_ok=True);stage=pathlib.Path(tempfile.mkdtemp(prefix=task.name+'.staging.',dir=task.parent));ACTIVE=stage
 payload={**expected,'status':'success','method':'TEMPL author core, same-row self-rerun','rmsd_A':float(checked['rmsd']),'success_le_2A':bool(checked['rmsd_within_threshold']),'posebusters_version':posebusters_version,'rdkit_version':rdBase.rdkitVersion,'metric':'PoseBusters check_rmsd, heavy_only=True, choose_by=rmsd','query_native_read_only_after_prediction_commit':True}
 atomic(stage/'result.json',payload);os.replace(stage,task);ACTIVE=None
 return {**expected,'status':'computed','rmsd_A':payload['rmsd_A'],'success_le_2A':payload['success_le_2A']}
def summarize():
 vals=[]
 for r in rows():
  p=OUT/f"{r['family']}__{r['query_id']}"/'result.json'
  if p.is_file():vals.append(json.loads(p.read_text()))
 import numpy as np
 d={'expected':44,'completed':len(vals),'success_rate':float(np.mean([v['success_le_2A'] for v in vals])) if vals else None,'mean_rmsd_A':float(np.mean([v['rmsd_A'] for v in vals])) if vals else None}
 atomic(OUT/'summary.json',d);return d
def main():
 p=argparse.ArgumentParser();p.add_argument('--index',type=int);p.add_argument('--list-count',action='store_true');p.add_argument('--summarize',action='store_true');p.add_argument('--fault-after-metric',action='store_true');p.add_argument('--hold-after-metric',type=float,default=0);a=p.parse_args()
 if a.list_count:print(44);return
 if a.summarize:print(json.dumps(summarize(),sort_keys=True));return
 if a.index is None:p.error('--index required')
 print(json.dumps(run(a.index,a.fault_after_metric,a.hold_after_metric),sort_keys=True))
if __name__=='__main__':main()
