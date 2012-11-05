import sys
import os
import unittest

sys.path.insert(0, os.path.dirname(__file__))
from pyteomics import biolccc

class TestPicklingFacilities(unittest.TestCase):
    def test_chemicalbasis_pickling(self):
        state = biolccc.rpAcnTfaChain.__getstate__()
        new_chembasis = biolccc.ChemicalBasis()
        new_chembasis.__setstate__(state)
        self.assertEqual(new_chembasis, biolccc.rpAcnTfaChain)

        import pickle
        unpickled_chembasis = pickle.loads(pickle.dumps(
            biolccc.ChemicalBasis(biolccc.RP_ACN_TFA_CHAIN)))
        self.assertEqual(unpickled_chembasis, biolccc.rpAcnTfaChain)

    def test_chromoconditions_pickling(self):
        state = biolccc.standardChromoConditions.__getstate__()
        new_chromatograph = biolccc.ChromoConditions()
        new_chromatograph.__setstate__(state)
        self.assertEqual(new_chromatograph, biolccc.standardChromoConditions)

        import pickle
        unpickled_chromoconditions = pickle.loads(pickle.dumps(
            biolccc.ChromoConditions()))
        self.assertEqual(unpickled_chromoconditions,
                         biolccc.standardChromoConditions)

if __name__ == '__main__':
    unittest.main()

