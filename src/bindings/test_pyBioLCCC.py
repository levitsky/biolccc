import sys
import os
import unittest

sys.path.insert(0, os.path.dirname(__file__))
import pyBioLCCC

class TestPicklingFacilities(unittest.TestCase):
    def test_chemicalbasis_pickling(self):
        state = pyBioLCCC.rpAcnTfaCoilBoltzmann.__getstate__()
        new_chembasis = pyBioLCCC.ChemicalBasis()
        new_chembasis.__setstate__(state)
        self.assertEqual(new_chembasis, pyBioLCCC.rpAcnTfaCoilBoltzmann)

    def test_chromoconditions_pickling(self):
        state = pyBioLCCC.standardChromoConditions.__getstate__()
        new_chromatograph = pyBioLCCC.ChromoConditions()
        new_chromatograph.__setstate__(state)
        self.assertEqual(new_chromatograph, pyBioLCCC.standardChromoConditions)

if __name__ == '__main__':
    unittest.main()

