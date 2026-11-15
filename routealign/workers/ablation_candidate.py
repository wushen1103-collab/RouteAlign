"""Controlled candidate stage for the identical-candidate ablation."""
import sys
from routealign.workers.fixed_pool_worker import main
if __name__=='__main__':
    sys.argv=[('--audit' if x=='--audit-only' else x) for x in sys.argv]
    sys.argv.insert(1,'candidate')
    main()
