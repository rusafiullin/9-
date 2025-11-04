import unittest
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main import Nuclear_Analyzer


class TestNuclearAnalyzer(unittest.TestCase):
    def setUp(self):
        self.analyzer = Nuclear_Analyzer()

    def test_validate_input_valid(self):
        Z, A = self.analyzer._validate_input(6, 12)
        self.assertEqual(Z, 6)
        self.assertEqual(A, 12)

    def test_validate_input_negative_Z(self):
        with self.assertRaises(ValueError):
            self.analyzer._validate_input(-1, 12)

    def test_validate_input_negative_A(self):
        with self.assertRaises(ValueError):
            self.analyzer._validate_input(6, -12)

    def test_validate_input_Z_greater_than_A(self):
        with self.assertRaises(ValueError):
            self.analyzer._validate_input(12, 6)

    def test_validate_input_invalid_types(self):
        with self.assertRaises(ValueError):
            self.analyzer._validate_input("abc", "def")

    def test_formula_known_values(self):
        energy = self.analyzer.formula(6, 12)
        self.assertIsInstance(energy, float)
        self.assertGreater(energy, 0)

    def test_formula_zero_A(self):
        energy = self.analyzer.formula(0, 0)
        self.assertEqual(energy, 0)

    def test_ud_energy_calculation(self):
        energy_per_nucleon = self.analyzer.ud_energy(26, 56)
        self.assertGreater(energy_per_nucleon, 0)
        self.assertLess(energy_per_nucleon, 10)

    def test_mass_calculation(self):
        mass = self.analyzer.mass(1, 1)
        self.assertIsInstance(mass, float)
        self.assertGreater(mass, 0)

    def test_radius_calculation(self):
        radius = self.analyzer.radius(92, 238)
        self.assertIsInstance(radius, float)
        self.assertGreater(radius, 0)

        radius_small = self.analyzer.radius(6, 12)
        radius_large = self.analyzer.radius(92, 238)
        self.assertGreater(radius_large, radius_small)

    def test_beta_stability(self):
        result = self.analyzer.beta(6, 12)
        self.assertIn(result, ["Устойчивое", "Неустойчивое бета-минус", "Неустойчивое бета-плюс"])

    def test_delenie_even_even(self):
        result = self.analyzer.delenie(92, 236)
        self.assertIn(result, ["Возможно", "Осколки нечетные"])

    def test_delenie_odd_A(self):
        result = self.analyzer.delenie(92, 235)
        self.assertEqual(result, "Нечетное A или Z")

    def test_edge_cases(self):
        mass_neutron_like = self.analyzer.mass(0, 1)
        self.assertIsInstance(mass_neutron_like, float)

        mass_proton = self.analyzer.mass(1, 1)
        self.assertIsInstance(mass_proton, float)

    def test_consistency(self):
        Z, A = 26, 56  # Железо-56

        energy = self.analyzer.formula(Z, A)
        energy_per_nucleon = self.analyzer.ud_energy(Z, A)
        mass = self.analyzer.mass(Z, A)

        self.assertAlmostEqual(energy_per_nucleon, energy / A, places=2)

        self.assertGreater(mass, 0)


class TestNuclearAnalyzerIntegration(unittest.TestCase):

    def setUp(self):
        self.analyzer = Nuclear_Analyzer()

    def test_complete_analysis(self):
        Z, A = 8, 16

        energy = self.analyzer.formula(Z, A)
        energy_per_nucleon = self.analyzer.ud_energy(Z, A)
        mass = self.analyzer.mass(Z, A)
        radius = self.analyzer.radius(Z, A)
        beta_result = self.analyzer.beta(Z, A)
        fission_result = self.analyzer.delenie(Z, A)

        self.assertIsInstance(energy, (int, float))
        self.assertIsInstance(energy_per_nucleon, (int, float))
        self.assertIsInstance(mass, (int, float))
        self.assertIsInstance(radius, (int, float))
        self.assertIsInstance(beta_result, str)
        self.assertIsInstance(fission_result, str)


if __name__ == '__main__':
    unittest.main(verbosity=2)