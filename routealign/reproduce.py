"""Recompute all reported strata and statistical comparisons from frozen observations."""
import os
os.environ.update(CUDA_VISIBLE_DEVICES='', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from routealign import settings


def compare(actual, expected, path='strata'):
    if isinstance(expected, dict):
        assert set(actual) == set(expected), path
        return sum(compare(actual[k], expected[k], path + '.' + k) for k in expected)
    if isinstance(expected, list):
        assert len(actual) == len(expected), path
        return sum(compare(a, b, path) for a, b in zip(actual, expected))
    if isinstance(expected, (float, int)) and not isinstance(expected, bool):
        assert math.isclose(actual, expected, rel_tol=1e-11, abs_tol=1e-12), (path, actual, expected)
        return 1
    assert actual == expected, (path, actual, expected)
    return 0


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=settings.workspace() / 'outputs/statistics')
    p.add_argument('--live', action='store_true', help='Read newly rerun worker outputs instead of frozen observations')
    args = p.parse_args()
    from routealign.workers import statistics, uniform_statistics
    ref = settings.PACKAGE_ROOT / 'reference'
    checks = json.loads((ref / 'checksums.json').read_text())
    for relative, digest in checks.items():
        assert hashlib.sha256((settings.PACKAGE_ROOT / relative).read_bytes()).hexdigest() == digest, relative
    raw = gzip.decompress((ref / 'observations.json.gz').read_bytes())
    data = json.loads(raw)
    args.output.mkdir(parents=True, exist_ok=True)
    results = []
    for label, module in [('hybrid', statistics)]:
        if not args.live:
            module.load = lambda label=label: (data[label + '_records'], hashlib.sha256(raw).hexdigest())
        result, _ = module.analyze()
        target = json.loads((ref / (label + '_summary.json')).read_text())
        compared = compare(result['strata'], target['strata'])
        (args.output / (label + '_summary.json')).write_text(json.dumps(result, indent=2) + '\n')
        rows = []
        for name, item in result['strata'].items():
            v = item['proposed']
            rows.append([name, item['n_pairs'], 'RouteAlign' if label == 'hybrid' else 'Uniform-medoid', v['rmsd_A_mean'], v['rmsd_A_std'], v['success_mean'], v['success_std'], v['PB20_mean'], v['PB20_std'], v['joint_mean'], v['joint_std'], 'self-rerun'])
            for method, v in item['baselines'].items():
                rows.append([name, item['n_pairs'], method, v['mean_rmsd_A'], v.get('rmsd_A_std'), v['success'], v.get('success_std'), v.get('PB20'), v.get('PB20_std'), v.get('joint_success'), v.get('joint_success_std'), 'published output' if method.startswith('Iktos') else 'self-rerun'])
        with (args.output / (label + '_table.csv')).open('w', newline='') as f:
            w = csv.writer(f); w.writerow(['stratum', 'pairs', 'method', 'mean_rmsd_A', 'std_rmsd_A', 'success', 'success_std', 'physical_pass', 'physical_pass_std', 'joint_success', 'joint_std', 'source']); w.writerows(rows)
        results.append({'analysis': label, 'strata': len(result['strata']), 'numeric_values_verified': compared})
    if not args.live:
        from routealign.ablation import analyze
        groups=json.loads(gzip.decompress((ref/'same_candidate_groups.json.gz').read_bytes()))
        ablation=analyze(groups)
        checked=compare(ablation,json.loads((ref/'ablation_summary.json').read_text()),'ablation')
        (args.output/'ablation_summary.json').write_text(json.dumps(ablation,indent=2)+'\n')
        results.append({'analysis':'same_candidate_ablation','strata':3,'numeric_values_verified':checked})
    print(json.dumps({'status': 'PASS', 'results': results, 'mode': 'live' if args.live else 'frozen'}, indent=2))


if __name__ == '__main__':
    main()
