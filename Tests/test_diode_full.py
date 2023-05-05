from unittest import TestCase
from unittest import TestCase
from diode import Diode
from sympy.physics.units import *
from constants import *
from pprint import pprint
import math
import sympy as sp


class TestDiode(TestCase):

    def setUp(self):
        self.diode = Diode()

        d = self.diode
        d.add_known_quantity('N_d', 1e17 / cm ** 3)
        d.add_known_quantity('N_a', 1e16 / cm ** 3)
        d.add_known_quantity('A', 2e-4 * cm ** 2)
        d.add_known_quantity('V_PN', 0.3 * V)

    def test_p_side(self):
        d = self.diode

        # Desired quantities
        # mu_nP = 1254.625
        # D_nP = 32.43476
        # tau_recP = 2.89775e-5
        # tau_genP = 2.17331e-3
        # L_n_rec_P = 3.0657450e-2
        # L_n_gen_P = 2.654935e-1

        pprint(d.get_known_quantities())
        self.assertAlmostEqual(remove_units(d.get_known_quantity('mu_nP')), 1254.625, delta=0.001 * 1254.625)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('D_nP')), 32.43476, delta=0.001 * 32.43476)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('tau_recP')), 2.89775e-5, delta=0.001 * 2.89775e-5)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('tau_genP')), 2.17331e-3, delta=0.001 * 2.17331e-3)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('L_n_rec_P')), 3.0657450e-2, delta=0.001 * 306.57450e-2)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('L_n_gen_P')), 2.654935e-1, delta=0.001 * 2.654935e-1)

    def test_n_side(self):
        d = self.diode

        # Desired quantities
        # mu_pN = 474.4022
        # D_pN = 12.26432
        # tau_recN = 1.253133e-5
        # tau_genN = 9.398496e-4
        # L_p_rec_N = 1.239710e-2
        # L_p_gen_N = 1.073589e-1

        pprint(d.get_known_quantities())
        self.assertAlmostEqual(remove_units(d.get_known_quantity('mu_pN')), 474.4022, delta=0.001 * 474.4022)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('D_pN')), 12.26432, delta=0.001 * 12.26432)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('tau_recN')), 1.253133e-5, delta=0.001 * 1.253133e-5)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('tau_genN')), 9.398496e-4, delta=0.001 * 9.398496e-4)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('L_p_rec_N')), 1.239710e-2, delta=0.001 * 1.239710e-2)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('L_p_gen_N')), 1.073589e-1, delta=0.001 * 1.073589e-1)

    def test_scr(self):
        d = self.diode

        # Desired quantities
        # tau_recSCR = 2.075443e-5
        # tau_genSCR = 1.556582e-3

        pprint(d.get_known_quantities())
        self.assertAlmostEqual(remove_units(d.get_known_quantity('tau_recSCR')), 2.075443e-5, delta=0.001 * 2.075443e-5)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('tau_genSCR')), 1.556582e-3, delta=0.001 * 1.556582e-3)

    def test_thermal_equilibrium(self):
        d = self.diode

        d.add_known_quantity('V_PN', 0 * volt)

        # Desired quantities
        # V_bi = 0.770350
        # W_dP = 3.009364e-5
        # W_dN = 3.009364e-6
        # W_d = 3.310300e-5
        # epsilon_x_max = -4.654260e4

        pprint(d.get_known_quantities())
        self.assertAlmostEqual(remove_units(d.get_known_quantity('V_bi')), 0.770350, delta=0.001 * 0.770350)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('W_dP')), 3.009364e-5, delta=0.001 * 3.009364e-5)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('W_dN')), 3.009364e-6, delta=0.001 * 3.009364e-6)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('W_d')), 3.310300e-5, delta=0.001 * 3.310300e-5)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('epsilon_x_max')), -4.654260e4, delta=0.001 * 4.654260e4)

    def test_forward_bias(self):
        d = self.diode
        d.add_known_quantity('V_PN', 0.3 * V)

        # Desired quantities
        # W_dP = 2.351480e-5
        # W_dN = 2.351480e-6
        # W_d = 2.586628e-5

        pprint(d.get_known_quantities())
        self.assertAlmostEqual(remove_units(d.get_known_quantity('W_d')), 2.586628e-5, delta=0.001 * 2.586628e-5)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('W_dP')), 2.351480e-5, delta=0.001 * 2.351480e-5)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('W_dN')), 2.351480e-6, delta=0.001 * 2.351480e-6)


    def test_saturation_currents(self):
        d = self.diode

        area = remove_units(d.get_known_quantity('A'))

        # Desired quantities
        # a * J_Sdiff = 3.881349e-16 + 3.629370e-17
        # a * J_Sscr = 2.13657e-13

        pprint(d.get_known_quantities())
        self.assertAlmostEqual(remove_units(d.get_known_quantity('J_Sdiff')) * area, 3.881349e-16 + 3.629370e-17, delta=0.001 * (3.881349e-16 + 3.629370e-17))
        self.assertAlmostEqual(remove_units(d.get_known_quantity('J_Sscr')) * area, 2.13657e-13, delta=0.001 * 2.13657e-13)

    def test_currents(self):
        d = self.diode

        a = remove_units(d.get_known_quantity('A'))

        # Desired quantities
        # a * J_nxdiffP = 4.253307e-11
        # a * J_pxdiffN = 3.977179e-12
        # I_Ddiff = 4.651025-11
        # a * J_Dscr = 7.051435e-11
        # I_D = 1.170246e-10

        pprint(d.get_known_quantities())
        self.assertAlmostEqual(remove_units(d.get_known_quantity('J_nxdiffP')) * a, 4.253307e-11, delta=0.001 * 4.253307e-11)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('J_pxdiffN')) * a, 3.977179e-12, delta=0.001 * 3.977179e-12)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('I_Ddiff')), 4.651025e-11, delta=0.001 * 4.651025e-11)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('J_Dscr')) * a, 7.051435e-11, delta=0.001 * 7.051435e-11)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('I_D')), 1.170246e-10, delta=0.001 * 1.170246e-10)

    def test_capacitances(self):
        d = self.diode

        d.add_known_quantity('V_PN', 0.3 * V)
        a = remove_units(d.get_known_quantity('A'))

        # Desired quantities
        # C_pn*dep = 8.009965e-12
        # C_pn*diffP * a = 4.767550e-14
        # C_pn*diffN * a = 1.927878e-15
        # C_pn*diff * a = 4.960338e-14
        # C_pn * a = 8.059568e-12

        pprint(d.get_known_quantities())
        self.assertAlmostEqual(remove_units(d.get_known_quantity('C_pn*dep_per_area')) * a, 8.009965e-12, delta=0.001 * 8.009965e-12)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('C_pn*diffP')) * a, 4.767550e-14, delta=0.001 * 4.767550e-14)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('C_pn*diffN')) * a, 1.927878e-15, delta=0.001 * 1.927878e-15)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('C_pn*diff')) * a, 4.960338e-14, delta=0.001 * 4.960338e-14)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('C_pn')) * a, 8.059568e-12, delta=0.001 * 8.059568e-12)

    def test_capacitances_equilibrium(self):
        d = self.diode

        d.add_known_quantity('V_PN', 0 * V)
        a = remove_units(d.get_known_quantity('A'))

        # Desired quantities
        # C_pn*dep = 6.258887e-12
        # C_pn*diff * a = 4.526508e-19

        pprint(d.get_known_quantities())
        self.assertAlmostEqual(remove_units(d.get_known_quantity('C_pn*dep_per_area')) * a, 6.258887e-12, delta=0.001 * 6.258887e-12)
        self.assertAlmostEqual(remove_units(d.get_known_quantity('C_pn*diff')) * a, 4.526508e-19, delta=0.001 * 4.526508e-19)

