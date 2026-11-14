import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from routealign import settings
from routealign.control import execute


class TaskControlTests(unittest.TestCase):
    def environment(self):
        return {**os.environ, 'CUDA_VISIBLE_DEVICES': '', 'OMP_NUM_THREADS': '1',
                'PYTHONPATH': str(settings.PACKAGE_ROOT)}

    @unittest.skipUnless(os.name == 'posix', 'Linux process-group control')
    def test_device_and_affinity_are_explicit(self):
        cpu = min(os.sched_getaffinity(0))
        code = "import os,json; print(json.dumps({'gpu':os.environ['CUDA_VISIBLE_DEVICES'],'cpu':sorted(os.sched_getaffinity(0))}))"
        rc, out, _ = execute([sys.executable, '-c', code], cwd=settings.PACKAGE_ROOT,
                              env=self.environment(), timeout=10, cpus=[cpu])
        self.assertEqual(rc, 0)
        self.assertEqual(json.loads(out), {'gpu': '', 'cpu': [cpu]})

    @unittest.skipUnless(os.name == 'posix', 'Linux process-group control')
    def test_timeout_stops_grandchild_before_commit(self):
        with tempfile.TemporaryDirectory() as d:
            target = Path(d) / 'must_not_exist'
            child = "import time,pathlib;time.sleep(3);pathlib.Path(" + repr(str(target)) + ").write_text('unexpected')"
            parent = "import subprocess,sys,time;subprocess.Popen([sys.executable,'-c'," + repr(child) + "]);time.sleep(20)"
            start = time.monotonic()
            rc, _, _ = execute([sys.executable, '-c', parent], cwd=d, env=self.environment(), timeout=.4)
            self.assertEqual(rc, 124)
            self.assertLess(time.monotonic() - start, 7)
            time.sleep(3.2)
            self.assertFalse(target.exists())


if __name__ == '__main__':
    unittest.main()
