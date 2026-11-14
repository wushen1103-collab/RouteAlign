"""Run one isolated task with a deadline and terminate its process group."""
import os
import signal
import subprocess
import sys
import threading

_ACTIVE = set()
_LOCK = threading.Lock()


def _terminate(proc, sig=signal.SIGTERM):
    if proc.poll() is not None:
        return
    try:
        if os.name == 'posix':
            os.killpg(proc.pid, sig)
        elif sig == signal.SIGTERM:
            proc.terminate()
        else:
            proc.kill()
    except ProcessLookupError:
        pass


def terminate_all():
    with _LOCK:
        processes = list(_ACTIVE)
    for proc in processes:
        _terminate(proc)


def execute(command, *, cwd, env, timeout, cpus=None):
    if cpus and os.name == 'posix':
        command = [sys.executable, '-m', 'routealign.task', ','.join(map(str, cpus)), *command]
    proc = subprocess.Popen(command, cwd=cwd, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True,
                            start_new_session=True)
    with _LOCK:
        _ACTIVE.add(proc)
    try:
        stdout, stderr = proc.communicate(timeout=timeout)
        return proc.returncode, stdout, stderr
    except (subprocess.TimeoutExpired, KeyboardInterrupt):
        _terminate(proc)
        try:
            stdout, stderr = proc.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            _terminate(proc, signal.SIGKILL)
            stdout, stderr = proc.communicate()
        return 124, stdout, stderr
    finally:
        with _LOCK:
            _ACTIVE.discard(proc)
