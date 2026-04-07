import contextlib
import io
import pathlib
import runpy
import unittest

from pyteomics import biolccc


EXAMPLES_DIR = pathlib.Path(__file__).resolve().parents[1] / 'src' / 'examples'

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


class TestDocumentedPythonExamples(unittest.TestCase):
    def run_example(self, filename):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            runpy.run_path(str(EXAMPLES_DIR / filename), run_name='__main__')
        return output.getvalue()

    def test_rt_calculation_example(self):
        output = self.run_example('rt_calculation.py')
        self.assertIn('The retention time of Ac-PEPTIDE-NH2 is', output)

    def test_chromoconditions_examples(self):
        output = self.run_example('chromoconditions.py')
        self.assertIn('The retention time of Ac-PEPTIDE-NH2 is', output)

        output = self.run_example('chromoconditions_dict.py')
        self.assertIn('columnLength', output)
        self.assertIn('The retention time of Ac-PEPTIDE-NH2 is', output)

    def test_mass_calculation_example(self):
        output = self.run_example('mass_calculation.py')
        self.assertIn('The average mass of Ac-PEPTIDE-NH2 is', output)
        self.assertIn('The monoisotopic mass of Ac-PEPTIDE-NH2 is', output)

    def test_kd_calculation_example(self):
        output = self.run_example('kd_calculation.py')
        self.assertIn('The coefficient of distribution of Ac-PEPTIDE-NH2 is', output)

    def test_gradient_example(self):
        output = self.run_example('gradient.py')
        self.assertIn('The retention time of Ac-PEPTIDE-NH2 in the custom gradient is', output)

    def test_advanced_chembasis_example(self):
        output = self.run_example('advanced_chembasis.py')
        self.assertIn('The bind energy of E is', output)
        self.assertIn('The bind energy of -NH2 is', output)
        self.assertIn('The retention time of Ac-PEhoPTIDE-NH2 is', output)
        self.assertIn('The monoisotopic mass of Ac-PEhoPTIDE-NH2 is', output)

    def test_sequence_parsing_example(self):
        output = self.run_example('sequence_parsing.py')
        lines = [line.strip() for line in output.splitlines() if line.strip()]
        self.assertEqual(lines[0], 'N-terminal hydrogen')
        self.assertIn('Glutamic acid', lines)
        self.assertEqual(lines[-1], 'C-terminal carboxyl group')

    def test_interpolation_example(self):
        output = self.run_example('interpolation.py')
        lines = [line.strip() for line in output.splitlines() if line.strip()]
        self.assertEqual(len(lines), 2)
        standard = float(lines[0])
        interpolated = float(lines[1])
        self.assertGreater(standard, 0.0)
        self.assertGreater(interpolated, 0.0)
        self.assertLess(abs(standard - interpolated), 0.5)


class TestCalculateRT(unittest.TestCase):
    def test_calculate_rt_with_default_conditions(self):
        rt = biolccc.calculateRT('PEPTIDE', biolccc.rpAcnTfaChain, biolccc.standardChromoConditions)
        self.assertAlmostEqual(rt, 11.0013889, places=5)

if __name__ == '__main__':
    unittest.main()

