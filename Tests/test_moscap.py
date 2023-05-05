from pprint import pprint
from unittest import TestCase
from sympy.physics.units import *
from moscap import MOSCAP
from constants import *

angstrom = 1e-8 * cm


class TestMOSCAP(TestCase):

    def setUp(self):
        self.moscap = MOSCAP()
        mc = self.moscap

        mc.add_known_quantity('N_aB', 5.0e17 * cm ** -3)
        mc.add_known_quantity('X_ox', 100 * angstrom)
        mc.add_known_quantity('Q_f', q * 5.0e10 * cm ** -2)
        mc.add_known_quantity('Q_it', 0 * cm ** -2)
        mc.add_known_quantity('N_d', 1.0e16 / cm ** 3)

        # TODO add area

    def test_parameters(self):
        mc = self.moscap

        # Desired quantities
        # N_a = 5.0e17
        # p_oP = 5.0e17
        # n_oP = 2.28980e2
        # X_ox = 100e-8
        # C_ox = 3.45313e-7
        # L_DB = 5.78196e-7
        # Q_f = 8.01088e-9
        # Q_it = 0
        # phi_PM = −0.96819
        # phi_FB = 0.45654

        pprint(mc.get_known_quantities())
        self.assertAlmostEqual(remove_units(mc.get_known_quantity('N_aB')), 5.0e17, delta=0.001 * 5.0e17)
        self.assertAlmostEqual(remove_units(mc.get_known_quantity('p_oP')), 5.0e17, delta=0.001 * 5.0e17)
        self.assertAlmostEqual(remove_units(mc.get_known_quantity('n_oP')), 2.28980e2, delta=0.001 * 2.28980e-2)
        self.assertAlmostEqual(remove_units(mc.get_known_quantity('X_ox')), 100e-8, delta=0.001 * 1e-8)
        self.assertAlmostEqual(remove_units(mc.get_known_quantity('C_ox')), 3.45313e-7, delta=0.001 * 3.45313e-7)
        self.assertAlmostEqual(remove_units(mc.get_known_quantity('L_DB')), 5.78196e-7, delta=0.001 * 5.78196e-7)
        self.assertAlmostEqual(remove_units(mc.get_known_quantity('Q_f')), 8.01088e-9, delta=0.001 * 8.01088e-9)
        self.assertAlmostEqual(remove_units(mc.get_known_quantity('Q_it')), 0, delta=0.001 * 0)
        self.assertAlmostEqual(remove_units(mc.get_known_quantity('phi_PM')), -0.96819, delta=0.001 * 0.96819)
        self.assertAlmostEqual(remove_units(mc.get_known_quantity('phi_FB')), 0.45654, delta=0.001 * 0.45654)

    def test_flatband_parameters(self):
        mc = self.moscap

        flatband = mc.get_known_quantity('V_FBN')
        mc.add_known_quantity('V_GB', flatband)

        # Desired quantities
        # epsilon_xOX@V_GB = −2.31989e4
        # V_ox = −0.023199
        # V_FBN = -0.991395
        # C_ox = 3.45313e-7

        pprint(mc.get_known_quantities())
        self.assertAlmostEqual(remove_units(mc.get_known_quantity('epsilon_xOX@V_GB')), -2.31989e4, delta=0.001 * 2.31989e4)
        self.assertAlmostEqual(remove_units(mc.get_known_quantity('V_ox')), -0.023199, delta=0.001 * 0.023199)
        self.assertAlmostEqual(remove_units(mc.get_known_quantity('V_FBN')), -0.991395, delta=0.001 * 0.991395)
        self.assertAlmostEqual(remove_units(mc.get_known_quantity('C_ox')), 3.45313e-7, delta=0.001 * 3.45313e-7)

    def test_threshold_parameters(self):
        mc = self.moscap

        self.fail()

    def test_with_VGB(self):
        mc = self.moscap

        mc.add_known_quantity('V_GB', 0.6 * V)

        self.fail()


        pprint(mc.get_known_quantities())

