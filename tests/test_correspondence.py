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

    def test_coordinate_stereochemistry_repair(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from routealign.stereo_repair import (
            infer_from_coordinates, isomeric_smiles, reflect_substituent,
            repair_candidate,
        )
        expected = Chem.AddHs(Chem.MolFromSmiles('C[C@H](O)C(=O)N'))
        self.assertEqual(AllChem.EmbedMolecule(expected, randomSeed=20261115), 0)
        expected = Chem.RemoveHs(expected)
        expected = infer_from_coordinates(expected)
        candidate = Chem.Mol(expected)
        centre = Chem.FindMolChiralCenters(
            candidate, includeUnassigned=False, includeCIP=True
        )[0][0]
        reflect_substituent(candidate, centre)
        self.assertNotEqual(
            isomeric_smiles(infer_from_coordinates(candidate)),
            isomeric_smiles(expected),
        )
        repaired, audit = repair_candidate(candidate, expected)
        self.assertEqual(isomeric_smiles(repaired), isomeric_smiles(expected))
        self.assertFalse(audit['native_pose_read'])
        self.assertGreater(audit['maximum_atom_movement_A'], 0.)

if __name__=='__main__':unittest.main()
