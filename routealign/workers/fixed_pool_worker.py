"""Isolated correction and same-candidate three-rule ablation.

Each task is an isolated, CPU-only subprocess with an outer deadline. Native
ligands are opened exclusively by the endpoint phase, after choices are saved.
"""
from pathlib import Path
import argparse,atexit,csv,hashlib,json,os,re,shutil,signal,subprocess,sys,tempfile,time
import numpy as np
from rdkit import Chem,rdBase
from routealign.atom_correspondence import matrix,self_test
import routealign.atom_correspondence as correspondence
from routealign import settings

R=settings.workspace()
OUT=R/'runs/same_candidate_ablation'
PREREG=R/'config/protocol.json'
MANIFEST=R/'config/cohort.csv'
DATA=R/'data/AlignDockBench'
SMINA=R/'tools/smina'
SELECTION_RULES=('medoid','receptor_consensus','routed')
ACTIVE=None

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def write(p,v):
    with Path(p).open('w') as f:json.dump(v,f,indent=2,sort_keys=True,allow_nan=False);f.flush();os.fsync(f.fileno())
def clean():
    global ACTIVE
    if ACTIVE is not None and ACTIVE.exists():shutil.rmtree(ACTIVE)
    ACTIVE=None
def stop(sig,_):clean();raise SystemExit(128+sig)
atexit.register(clean);signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
def checked(p):
    p.resolve().relative_to(R.resolve())
    if not p.is_file():raise FileNotFoundError(p)
    return p
def specs():
    p=load(PREREG)
    rows={f"{x['family']}__{x['query_id']}":x for x in csv.DictReader(MANIFEST.open()) if x['final_primary_eligible']=='True'}
    assert len(p['pairs'])==44 and len(p['outer_seed_groups'])==5
    return p,rows
def single(p):
    ms=list(Chem.SDMolSupplier(str(p),removeHs=True,sanitize=True))
    if len(ms)!=1 or ms[0] is None:raise ValueError('invalid single SDF '+str(p))
    return ms[0]
def source(pair,seed):
    low=pair['hard_low_2d'];name=pair['pair']
    repair=R/'runs'/('bounded_repair' if low else 'receptor_repair')/name/f'seed_{seed}'
    gen=R/'runs'/('graft' if low else 'receptor_score')/name/f'seed_{seed}'
    return repair/'selected.sdf',repair/('result.json' if low else 'metadata.json'),gen/'metadata.json',gen/('pose_pocket_hybrid.sdf' if low else 'selected.sdf')
def identity(phase,index,inputs):
    return {'phase':phase,'index':index,'worker_sha256':sha(__file__),
      'mapping_sha256':sha(Path(correspondence.__file__)),
      'preregister_sha256':sha(PREREG),'manifest_sha256':sha(MANIFEST),
      'rdkit_version':rdBase.rdkitVersion,'inputs':{str(checked(p).relative_to(R)):sha(p) for p in inputs}}
def begin(phase,index,expected,args):
    global ACTIVE
    dest=OUT/phase/f'{index:04d}'
    if (dest/'result.json').is_file():
        old=load(dest/'result.json')
        if old['identity']!=expected:raise RuntimeError('resume identity drift')
        if old['status']!='complete':raise RuntimeError('incomplete result')
        for f,h in old.get('outputs',{}).items():
            if sha(dest/f)!=h:raise RuntimeError('resume output digest drift')
        return dest,None,old
    if args.audit:return dest,None,{'status':'audit_pass','identity':expected,'cpu_affinity':sorted(os.sched_getaffinity(0)),'gpu_hidden':os.environ.get('CUDA_VISIBLE_DEVICES')==''}
    dest.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix=dest.name+'.staging.',dir=dest.parent));ACTIVE=stage
    return dest,stage,None
def finish(dest,stage,expected,payload,args):
    global ACTIVE
    if args.hold:time.sleep(args.hold)
    if args.fault:raise SystemExit(97)
    payload={'identity':expected,'status':'complete',**payload}
    write(stage/'result.json',payload);os.replace(stage,dest);ACTIVE=None
    return payload

def candidate(index,args):
    p,rows=specs();seeds=sum(p['outer_seed_groups'],[]);pair=p['pairs'][index//25];seed=seeds[index%25];row=rows[pair['pair']]
    pose,rep,gen,prerepair=source(pair,seed);protein=DATA/row['family']/row['template_id']/'protein.pdb'
    expected=identity('candidate',index,[pose,rep,gen,prerepair,protein,SMINA]);dest,stage,old=begin('candidate',index,expected,args)
    if old:return old
    g=load(gen);repair=load(rep)
    if pair['hard_low_2d']:
        cmd=[str(SMINA),'-r',str(protein),'-l',str(prerepair),'--score_only','--cpu','1','--seed',str(seed)]
        done=subprocess.run(cmd,cwd=stage,capture_output=True,text=True,timeout=120)
        hits=re.findall(r'Affinity:\s*([-+0-9.eE]+)',done.stdout)
        if done.returncode or len(hits)!=1:raise RuntimeError('smina score_only failed '+done.stdout[-1000:])
        selected=g['selected']['pocket_hybrid'];h=float(selected['shape_distance'])
        strain=(float(selected['hybrid_score'])-h)/.02
        if strain < -1e-7:raise RuntimeError('negative relative strain')
        receptor=float(hits[0])+.10*max(0,strain)+.05*h
        provenance={'affinity':float(hits[0]),'relative_strain':strain,'shape_distance':h,'new_score_only':True}
        archived_pb=bool(repair['PB20_pass'])
    else:
        selected=[x for x in g['scores'] if x['unique_index']==g['selected_unique_index']]
        assert len(selected)==1
        receptor=float(selected[0]['selection_score']);provenance={'new_score_only':False,'saved_selected_score':selected[0]}
        archived_pb=bool(repair['final_PB_pass'])
    # Re-run every candidate's physical checks once, sharing them across rules.
    from posebusters import PoseBusters,__version__
    np.random.seed(20261115);from rdkit import rdBase as rb;rb.SeedRandomNumberGenerator(20261115)
    frame=PoseBusters(config='dock',max_workers=0).bust(str(pose),mol_cond=str(protein),full_report=False)
    assert len(frame)==1
    checks={str(k):bool(v) for k,v in frame.iloc[0].to_dict().items()}
    science={k:v for k,v in checks.items() if k not in ['mol_pred_loaded','mol_cond_loaded']};assert len(science)==19
    if not checks.get('mol_pred_loaded') or not checks.get('mol_cond_loaded'):raise RuntimeError('PB load failed')
    pb=all(science.values())
    return finish(dest,stage,expected,{'pair':pair['pair'],'seed':seed,'pose':str(pose.relative_to(R)),
        'pose_sha256':sha(pose),'prerepair_sha256':sha(prerepair),'receptor_score':receptor,'score_provenance':provenance,
        'PB_pass':pb,'checks':science,'archived_PB_pass':archived_pb,'posebusters_version':__version__,
        'native_read':False,'budget_shared_across_rules':True},args)

def select(index,args):
    p,rows=specs();pair=p['pairs'][index//5];gi=index%5;seeds=p['outer_seed_groups'][gi]
    inputs=[OUT/'candidate'/f'{index//5*25+gi*5+i:04d}'/'result.json' for i in range(5)]
    cs=[load(x) for x in inputs];poses=[R/x['pose'] for x in cs]
    expected=identity('select',index,inputs+poses);dest,stage,old=begin('select',index,expected,args)
    if old:return old
    mols=[single(f) for f in poses];d,mapping=matrix(mols)
    # Explicit reverse-direction agreement guards graph matcher asymmetry.
    from routealign.atom_correspondence import distance
    for i in range(5):
        for j in range(i+1,5):
            reverse,_=distance(mols[j],mols[i]);assert abs(reverse-d[i,j])<1e-9
    med=np.array([np.median(np.delete(d[i],i)) for i in range(5)])
    def z(a):
        a=np.asarray(a,float);return np.zeros_like(a) if a.std()<=1e-12 else (a-a.mean())/a.std()
    objective=z([c['receptor_score'] for c in cs])+5*z(med)
    eligible=[i for i,c in enumerate(cs) if c['PB_pass']] or list(range(5))
    mi=min(range(5),key=lambda i:(med[i],float(np.mean(np.delete(d[i],i))),seeds[i]));ri=min(eligible,key=lambda i:(objective[i],seeds[i]))
    selected={'medoid':mi,'receptor_consensus':ri,'routed':mi if pair['hard_low_2d'] else ri}
    arrays=[np.asarray(m.GetConformer().GetPositions()) for m in mols]
    oldd=np.array([[np.sqrt(np.mean(np.sum((x-y)**2,axis=1))) for y in arrays] for x in arrays]);oldmed=[np.median(np.delete(oldd[i],i)) for i in range(5)]
    oldobj=z([c['receptor_score'] for c in cs])+5*z(oldmed)
    oldmi=min(range(5),key=lambda i:(oldmed[i],float(np.mean(np.delete(oldd[i],i))),seeds[i]));oldri=min(eligible,key=lambda i:(oldobj[i],seeds[i]))
    outputs={};choices={}
    for rule,i in selected.items():
        filename=rule+'.sdf';shutil.copy2(poses[i],stage/filename);outputs[filename]=sha(stage/filename)
        choices[rule]={'index':i,'seed':seeds[i],'pose_sha256':cs[i]['pose_sha256'],'candidate_record':str(inputs[i].relative_to(R)),
                       'legacy_index_seed':seeds[oldmi if rule=='medoid' or (rule=='routed' and pair['hard_low_2d']) else oldri]}
    return finish(dest,stage,expected,{'pair':pair['pair'],'group':gi+1,'hard_low_2d':pair['hard_low_2d'],
        'hard_flex':pair['hard_flex'],'hard_scaffold':pair['hard_scaffold'],'candidate_seeds':seeds,
        'candidate_hashes':[c['pose_sha256'] for c in cs],'distance_matrix':d.tolist(),'mappings':mapping,
        'medians':med.tolist(),'receptor_scores':[c['receptor_score'] for c in cs],'objective':objective.tolist(),
        'PB_pass':[c['PB_pass'] for c in cs],'eligible':eligible,'choices':choices,'outputs':outputs,
        'native_read':False,'same_candidate_pool_for_all_rules':True,'mapping_fitting':False,'medoid_ties':'median, mean, seed'},args)

def endpoint(index,args):
    p,rows=specs();pair=p['pairs'][index//5];row=rows[pair['pair']];native=DATA/row['family']/row['query_id']/'ligand.sdf'
    selection=OUT/'select'/f'{index:04d}'/'result.json';sel=load(selection)
    inputs=[selection,native]+[selection.parent/(r+'.sdf') for r in SELECTION_RULES]
    inputs += [R/sel['choices'][r]['candidate_record'] for r in SELECTION_RULES]
    expected=identity('endpoint',index,inputs);dest,stage,old=begin('endpoint',index,expected,args)
    if old:return old
    from posebusters.modules.rmsd import check_rmsd
    results={}
    for rule in SELECTION_RULES:
        pose=selection.parent/(rule+'.sdf');assert sha(pose)==sel['outputs'][rule+'.sdf']
        metric=check_rmsd(single(pose),single(native),rmsd_threshold=2.,heavy_only=True,choose_by='rmsd')['results']
        rmsd=float(metric['rmsd']);success=bool(metric['rmsd_within_threshold']);c=load(R/sel['choices'][rule]['candidate_record'])
        results[rule]={'rmsd_A':rmsd,'success':success,'PB_pass':c['PB_pass'],'joint':success and c['PB_pass'],
                      'seed':sel['choices'][rule]['seed'],'checks':c['checks'],'pose_sha256':sha(pose)}
    return finish(dest,stage,expected,{'pair':sel['pair'],'group':sel['group'],'hard_low_2d':sel['hard_low_2d'],
        'hard_flex':sel['hard_flex'],'hard_scaffold':sel['hard_scaffold'],'results':results,
        'native_access':'post-selection endpoint only','selection_sha256':sha(selection)},args)

def batch(phase,workers):
    from concurrent.futures import ThreadPoolExecutor,as_completed
    n=1100 if phase=='candidate' else 220;cpus=sorted(os.sched_getaffinity(0))
    if workers>min(8,len(cpus)-30):raise RuntimeError('CPU reserve/concurrency exceeded')
    if not (OUT/'control_audit.json').is_file():raise RuntimeError('control audit required before batching')
    def task(i):
        cmd=['taskset','-c',str(cpus[i%workers]),'timeout','-k','5s',('170s' if phase=='candidate' else '55s'),sys.executable,'-B',__file__,phase,'--index',str(i)]
        run=subprocess.run(cmd,start_new_session=True,capture_output=True,timeout=185 if phase=='candidate' else 70)
        (OUT/'logs').mkdir(exist_ok=True)
        (OUT/'logs'/f'{phase}_{i:04d}.log').write_bytes(run.stdout+run.stderr)
        if run.returncode:raise RuntimeError(f'{phase} {i}: exit {run.returncode} '+run.stderr.decode(errors='replace')[-1200:])
        return i
    t=time.time();completed=0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for future in as_completed([pool.submit(task,i) for i in range(n)]):
            future.result();completed+=1
            if completed%25==0 or completed==n:print(json.dumps({'phase':phase,'completed':completed,'total':n,'elapsed_s':round(time.time()-t,1)}),flush=True)
    write(OUT/(phase+'_batch_receipt.json'),{'n':n,'workers':workers,'reserve_min':30,'elapsed_s':time.time()-t})

def main():
    p=argparse.ArgumentParser();p.add_argument('phase',choices=['self-test','candidate','select','endpoint','batch']);p.add_argument('--index',type=int,default=0);p.add_argument('--audit',action='store_true');p.add_argument('--fault',action='store_true');p.add_argument('--hold',type=float,default=0);p.add_argument('--kind',choices=['candidate','select','endpoint']);p.add_argument('--workers',type=int,default=4);a=p.parse_args()
    assert Path.cwd().resolve()==R.resolve() and os.environ.get('CUDA_VISIBLE_DEVICES')==''
    OUT.resolve().relative_to(R.resolve());OUT.mkdir(exist_ok=True)
    if a.phase=='self-test':print(json.dumps(self_test()));return
    if a.phase=='batch':batch(a.kind,a.workers);return
    if not a.audit and os.getsid(0)!=os.getpgrp():raise RuntimeError('task requires isolated session/process group')
    result={'candidate':candidate,'select':select,'endpoint':endpoint}[a.phase](a.index,a)
    print(json.dumps({'status':result['status'],'phase':a.phase,'index':a.index,'choices':result.get('choices'),'results':result.get('results'),'PB_pass':result.get('PB_pass')},sort_keys=True))
if __name__=='__main__':main()
