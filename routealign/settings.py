"""Locations and resource choices supplied by the operator."""
from pathlib import Path
import os
import shutil

PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def workspace():
    return Path(os.environ.get('ROUTEALIGN_WORKSPACE', PACKAGE_ROOT)).resolve()


def cpu_index():
    available = sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else [0]
    chosen = int(os.environ.get('ROUTEALIGN_CPU', available[0]))
    if chosen not in available:
        raise RuntimeError('Requested CPU is outside the allowed affinity mask')
    return chosen


def openbabel():
    binary = os.environ.get('ROUTEALIGN_OBABEL') or shutil.which('obabel')
    if not binary:
        raise RuntimeError('Install Open Babel or set ROUTEALIGN_OBABEL to its executable')
    return Path(binary).resolve()
