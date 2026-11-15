"""Recompute the same-candidate selector ablation from frozen observations.

This is a statistics-level replay. Candidate graph distances and physical
checks were computed before selection; native-pose endpoints are accessed
only after the selected candidate index has been fixed.
"""
import os
os.environ.update(CUDA_VISIBLE_DEVICES='',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
import argparse,gzip,json
from pathlib import Path
import numpy as np
from routealign import settings

RULES=('medoid','receptor_consensus','routed')

def choose(features):
    d=np.asarray(features['distance_matrix'],float)
    if d.shape!=(5,5) or not np.allclose(d,d.T,atol=1e-12):raise ValueError('invalid five-candidate distances')
    seeds=features['seeds']
    med=np.asarray([np.median(np.delete(d[i],i)) for i in range(5)])
    avg=np.asarray([np.mean(np.delete(d[i],i)) for i in range(5)])
    def z(v):
        v=np.asarray(v,float);return np.zeros_like(v) if v.std()<=1e-12 else (v-v.mean())/v.std()
    objective=z(features['receptor_scores'])+5*z(med)
    eligible=[i for i,p in enumerate(features['physical_pass']) if p] or list(range(5))
    mi=min(range(5),key=lambda i:(med[i],avg[i],seeds[i]))
    ri=min(eligible,key=lambda i:(objective[i],seeds[i]))
    return {'medoid':mi,'receptor_consensus':ri,'routed':mi if features['low_similarity'] else ri}

def analyze(groups):
    from routealign.workers.statistics import exact_family_signflip,bootstrap,holm_adjust
    pairs={}
    for g in groups:
        choices=choose(g['features'])
        row=pairs.setdefault(g['pair'],{'pair':g['pair'],'family':g['family'],'low':g['features']['low_similarity'],
            'rules':{r:{'rmsd':[],'pb':[],'joint':[]} for r in RULES}})
        for rule,i in choices.items():
            # Endpoint evidence is intentionally read after choices are fixed.
            rmsd=g['endpoints']['rmsd_A'][i];pb=g['features']['physical_pass'][i]
            r=row['rules'][rule];r['rmsd'].append(rmsd);r['pb'].append(pb);r['joint'].append(rmsd<=2 and pb)
    if len(pairs)!=44 or any(len(r['rules']['routed']['rmsd'])!=5 for r in pairs.values()):raise ValueError('incomplete frozen cohort')
    result={};raw=[]
    for name,rows in [('all',list(pairs.values())),('low',[r for r in pairs.values() if r['low']]),('other',[r for r in pairs.values() if not r['low']])]:
        entry={'n_pairs':len(rows),'n_families':len(set(r['family'] for r in rows)),'rules':{},'comparisons':{}}
        for rule in RULES:
            rmsd=np.asarray([r['rules'][rule]['rmsd'] for r in rows]);metrics={}
            for metric,a in [('rmsd_A',rmsd),('success',rmsd<=2),('PB_pass',np.asarray([r['rules'][rule]['pb'] for r in rows])),('joint',np.asarray([r['rules'][rule]['joint'] for r in rows]))]:
                v=a.mean(0);metrics[metric]={'mean':float(v.mean()),'std':float(v.std(ddof=1)),'per_group':v.tolist()}
            entry['rules'][rule]=metrics
        proposed=np.asarray([r['rules']['routed']['rmsd'] for r in rows])
        for rule in RULES[:2]:
            baseline=np.asarray([r['rules'][rule]['rmsd'] for r in rows]);delta=proposed.mean(1)-baseline.mean(1)
            test=exact_family_signflip(rows,delta);raw.append((name,rule,test['one_sided_p']))
            entry['comparisons'][rule]={'delta_mean_A':float(delta.mean()),'family_signflip':test,
                'hierarchical_bootstrap':bootstrap(rows,proposed,baseline),'pair_wins':int((delta < -1e-12).sum()),
                'pair_ties':int((abs(delta)<=1e-12).sum()),'pair_losses':int((delta>1e-12).sum())}
        result[name]=entry
    adjusted=holm_adjust(raw)
    for name,rule,_ in raw:result[name]['comparisons'][rule]['family_signflip']['holm_across_6']=adjusted[name,rule]
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=settings.workspace()/'outputs/statistics/ablation_summary.json');a=p.parse_args()
    ref=settings.PACKAGE_ROOT/'reference'
    groups=json.loads(gzip.decompress((ref/'same_candidate_groups.json.gz').read_bytes()))
    result=analyze(groups)
    from routealign.reproduce import compare
    n=compare(result,json.loads((ref/'ablation_summary.json').read_text()),'ablation')
    a.output.parent.mkdir(exist_ok=True,parents=True);a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','groups':len(groups),'candidates':1100,'numeric_values_verified':n,'output':str(a.output)}))

if __name__=='__main__':main()
