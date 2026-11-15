"""Portable entry point for archived scientific workers (Linux)."""
import argparse
import csv
import json
import os
from pathlib import Path
import sys
import signal
import threading
import queue
from importlib.metadata import version
from concurrent.futures import ThreadPoolExecutor
from routealign import settings
from routealign.control import execute, terminate_all

STAGES = {
    'ablation_candidate':1100, 'ablation_select':220, 'ablation_evaluate':220,
    'graft': 1100, 'repair': 1100, 'bounded_repair': 1100,
    'medoid_select': 220, 'medoid_evaluate': 220,
    'receptor_score': 450, 'receptor_repair': 450,
    'receptor_select': 90, 'receptor_evaluate': 90,
    'templ_generate': 44, 'templ_evaluate': 44,
    'lsalign_generate': 1100, 'lsalign_evaluate': 1100, 'lsalign_validate': 220,
    'roshambo_generate': 220, 'roshambo_evaluate': 220, 'roshambo_validate': 220,
}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage', required=True, choices=STAGES)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--index', type=int)
    g.add_argument('--indices', help='Comma-separated task indices')
    g.add_argument('--all', action='store_true')
    p.add_argument('--workspace', type=Path, default=settings.workspace())
    p.add_argument('--workers', type=int, default=1)
    p.add_argument('--timeout', type=float, default=900)
    p.add_argument('--gpu', type=int, help='Explicit physical device for ROSHAMBO generation only')
    p.add_argument('--audit-only', action='store_true')
    args, extra = p.parse_known_args()
    if os.name != 'posix':
        p.error('Scientific workers require Linux; frozen statistics also run on other platforms')
    required_rdkit = '2024.9.6' if args.stage == 'templ_generate' else '2023.9.6' if args.stage == 'roshambo_generate' else '2026.9.1'
    if version('rdkit') != required_rdkit:
        p.error(f'{args.stage} requires RDKit {required_rdkit}; use the matching generation or endpoint environment')
    if args.stage in ['ablation_candidate', 'ablation_evaluate', 'repair', 'receptor_repair', 'medoid_evaluate', 'receptor_evaluate', 'templ_evaluate', 'lsalign_evaluate', 'lsalign_validate', 'roshambo_evaluate', 'roshambo_validate'] and version('posebusters') != '0.5.1':
        p.error('The reported endpoint protocol requires PoseBusters 0.5.1')
    cpus = sorted(os.sched_getaffinity(0))
    cap = max(1, len(cpus) - 30)
    if not 1 <= args.workers <= cap:
        p.error(f'Use between 1 and {cap} workers to retain the CPU reserve')
    if args.timeout <= 0:
        p.error('Timeout must be positive')
    if args.stage == 'roshambo_generate' and (args.gpu is None or args.workers != 1):
        p.error('ROSHAMBO requires an explicit --gpu and one worker')
    elif args.stage != 'roshambo_generate' and args.gpu is not None:
        p.error('--gpu is only applicable to ROSHAMBO generation')
    indices = list(range(STAGES[args.stage])) if args.all else [args.index] if args.index is not None else [int(x) for x in args.indices.split(',')]
    if not indices or any(not 0 <= i < STAGES[args.stage] for i in indices):
        p.error('Task index outside the stage range')
    if len(set(indices)) != len(indices):
        p.error('Duplicate task indices are not allowed')
    root = args.workspace.resolve()
    if not (root / 'config/protocol.json').is_file():
        p.error('Workspace lacks config/protocol.json; run routealign.prepare first')
    env = dict(os.environ)
    env.update(ROUTEALIGN_WORKSPACE=str(root), CUDA_VISIBLE_DEVICES='' if args.gpu is None else str(args.gpu),
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               NUMEXPR_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    env['PYTHONPATH'] = str(settings.PACKAGE_ROOT) + os.pathsep + env.get('PYTHONPATH', '')
    if args.gpu is not None:
        env['ROUTEALIGN_GPU'] = str(args.gpu)
    rows = list(csv.DictReader((root / 'config/cohort.csv').open()))
    rows = sorted((r for r in rows if r['final_primary_eligible'] == 'True'), key=lambda r: (r['family'], r['query_id']))
    slots = queue.Queue()
    for cpu in cpus[:args.workers]: slots.put(cpu)
    aborted = threading.Event()
    def stop(signum, _frame):
        aborted.set()
        terminate_all()
        raise SystemExit(128 + signum)
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    def task(item):
        position, index = item
        if aborted.is_set(): return 125
        cpu = slots.get()
        if aborted.is_set():
            slots.put(cpu)
            return 125
        local_env = dict(env, ROUTEALIGN_CPU=str(cpu))
        worker = args.stage
        if worker == 'graft' and (rows[index // 25]['family'], rows[index // 25]['query_id']) == ('2ica_2IC_A_1', '3m6f'):
            worker = 'graft_fallback'
        cmd = [sys.executable, '-B', '-m', 'routealign.workers.' + worker, '--index', str(index)]
        if args.audit_only:
            if worker in ['graft', 'graft_fallback', 'templ_generate']:
                cmd = cmd[:-2] + ['--list']
            elif worker == 'templ_evaluate':
                cmd = cmd[:-2] + ['--list-count']
            else:
                cmd.append('--audit-only')
        cmd.extend(extra)
        try:
            result = execute(cmd, cwd=root, env=local_env, timeout=args.timeout, cpus=[cpu])
        finally:
            slots.put(cpu)
        if result[0]:
            aborted.set()
            terminate_all()
        log = root / 'logs' / args.stage
        log.mkdir(parents=True, exist_ok=True)
        (log / f'{index}.out').write_text(result[1], encoding='utf-8')
        (log / f'{index}.err').write_text(result[2], encoding='utf-8')
        print(json.dumps({'stage': args.stage, 'index': index, 'exit_code': result[0]}), flush=True)
        return result[0]
    # Separate worker processes each receive a single CPU; the executor does no chemistry.
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        codes = list(pool.map(task, enumerate(indices)))
    raise SystemExit(0 if all(code == 0 for code in codes) else next(code for code in codes if code))


if __name__ == '__main__':
    main()
