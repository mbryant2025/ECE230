from unittest import TestCase
from diode import Diode
from sympy.physics.units import *
from constants import *
from pprint import pprint


class TestDiode(TestCase):

    def setUp(self):
        self.diode = Diode()

    def check_Nd_params(self):
        d = self.diode
        # Desired values
        # N_d = 1e16
        # mu_nN = 1247.98792
        # mu_pN = 498.45997
        # D_pN = 12.88626
        # tau_recN = 1.27910e-4
        # tau_genN = 9.59325e-3
        # L_p_rec_N = 405.99033 * 1e-4
        # L_p_gen_N = 3515.97938 * 1e-4
        pprint(d.get_known_quantities())
        self.assertAlmostEqual(remove_units(d.get_known_quantity('N_d')), 1e16, delta=0.001*1e16)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('mu_nN')), 1247.98792, delta=0.001*1247.98792)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('mu_pN')), 498.45997, delta=0.001*498.45997)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('D_pN')), 12.88626, delta=0.001*12.88626)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('tau_recN')), 1.27910e-4, delta=0.001*1.27910e-4)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('tau_genN')), 9.59325e-3, delta=0.001*9.59325e-3)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('L_p_rec_N')), 405.99033 * 1e-4, delta=0.001*405.99033 * 1e-4)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('L_p_gen_N')), 3515.97938 * 1e-4, delta=0.001*3515.97938 * 1e-4)

    def check_Na_params(self):
        d = self.diode
        # Desired values
        # N_a = 1e16
        # mu_nP = 1254.62510
        # mu_pP = 415.47881
        # D_nP = 32.43476
        # tau_recP = 2.89775 e-5
        # tau_genP = 2.17331 e-3
        # L_n_rec_P = 306.57450 * 1e-4
        # L_n_gen_P = 2655.01305 * 1e-4
        pprint(d.get_known_quantities())
        self.assertAlmostEqual(remove_units(d.get_known_quantity('N_a')), 1e16, delta=0.001*1e16)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('mu_nP')), 1254.62510, delta=0.001*1254.62510)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('mu_pP')), 415.47881, delta=0.001*415.47881)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('D_nP')), 32.43476, delta=0.001*32.43476)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('tau_recP')), 2.89775e-5, delta=0.001*2.89775e-5)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('tau_genP')), 2.17331e-3, delta=0.001*2.17331e-3)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('L_n_rec_P')), 306.57450 * 1e-4, delta=0.001*306.57450 * 1e-4)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('L_n_gen_P')), 2655.01305 * 1e-4, delta=0.001*2655.01305 * 1e-4)

    def test_given_N_d(self):
        d = self.diode
        d.add_known_quantity('N_d', 1e16 / cm**3)
        self.check_Nd_params()

    def test_given_mu_nN(self):
        d = self.diode
        d.add_known_quantity('mu_nN', 1247.98792 * cm ** 2 / (V * s))
        # add everything but N_d
        d.add_known_quantity('N_d', None)
        d.add_known_quantity('D_pN', 12.88626 * cm**2 / s)
        d.add_known_quantity('tau_recN', 1.27910e-4 * s)
        d.add_known_quantity('tau_genN', 9.59325e-3 * s)
        d.add_known_quantity('L_p_rec_N', 405.99033 * 1e-4 * cm)
        d.add_known_quantity('L_p_gen_N', 3515.97938 * 1e-4 * cm)
        self.check_Nd_params()

    def test_given_mu_pN(self):
        d = self.diode
        d.add_known_quantity('mu_pN', 498.45997 * cm**2 / (V * s))
        # add everything but N_d
        d.add_known_quantity('N_d', None)
        d.add_known_quantity('D_pN', 12.88626 * cm**2 / s)
        d.add_known_quantity('tau_recN', 1.27910e-4 * s)
        d.add_known_quantity('tau_genN', 9.59325e-3 * s)
        d.add_known_quantity('L_p_rec_N', 405.99033 * 1e-4 * cm)
        d.add_known_quantity('L_p_gen_N', 3515.97938 * 1e-4 * cm)
        self.check_Nd_params()

    def test_given_N_a(self):
        d = self.diode
        d.add_known_quantity('N_a', 1e16 / cm**3)
        self.check_Na_params()

    def test_given_mu_nP(self):
        d = self.diode
        d.add_known_quantity('mu_nP', 1254.62510 * cm ** 2 / (V * s))
        # add everything but N_a
        d.add_known_quantity('N_a', None)
        d.add_known_quantity('D_nP', 32.43476 * cm**2 / s)
        d.add_known_quantity('tau_recP', 2.89775e-5 * s)
        d.add_known_quantity('tau_genP', 2.17331e-3 * s)
        d.add_known_quantity('L_n_rec_P', 306.57450 * 1e-4 * cm)
        d.add_known_quantity('L_n_gen_P', 2655.01305 * 1e-4 * cm)
        self.check_Na_params()

    def test_given_mu_pP(self):
        d = self.diode
        d.add_known_quantity('mu_pP', 415.47881 * cm**2 / (V * s))
        # add everything but N_a
        d.add_known_quantity('N_a', None)
        d.add_known_quantity('D_nP', 32.43476 * cm**2 / s)
        d.add_known_quantity('tau_recP', 2.89775e-5 * s)
        d.add_known_quantity('tau_genP', 2.17331e-3 * s)
        d.add_known_quantity('L_n_rec_P', 306.57450 * 1e-4 * cm)
        d.add_known_quantity('L_n_gen_P', 2655.01305 * 1e-4 * cm)
        self.check_Na_params()





