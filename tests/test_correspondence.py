"""Scientific invariants for candidate atom correspondence."""
import unittest
try:
    import rdkit
except ImportError:
    rdkit=None

@unittest.skipIf(rdkit is None,'install chemistry extras for graph-mapping tests')
class CorrespondenceTests(unittest.TestCase):
    def test_symmetry_renumbering_and_receptor_frame(self):
        from routealign.atom_correspondence import self_test
        checks=self_test()
        self.assertEqual(checks['resonance_terminal_swap']['rmsd'],0.)
        self.assertEqual(checks['c1ccccc1']['maps'],12)
        self.assertEqual(checks['CC(=O)[O-]']['translated_rmsd'],3.)
        self.assertTrue(checks['stereo_tag_drift_recorded_without_coordinate_change'])

if __name__=='__main__':unittest.main()
