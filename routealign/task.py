"""Set task affinity in the child, before importing scientific libraries."""
import os
import sys

if __name__ == '__main__':
    cpus = {int(x) for x in sys.argv[1].split(',')}
    if not cpus <= os.sched_getaffinity(0):
        raise SystemExit('Task affinity exceeds the allowed CPU mask')
    os.sched_setaffinity(0, cpus)
    os.execvpe(sys.argv[2], sys.argv[2:], os.environ)
