from ece230 import ECE230
from constants import *
import sympy as sp
from sympy.physics.units import *


class Diode(ECE230):

    def __init__(self, T = 300 * kelvin):
        super().__init__(T)
        self.calculation_functions.extend([self.calc_V_bi, self.calc_V_bin, self.calc_V_bip, self.calc_W_dN, self.calc_W_dP,
                                            self.calc_po_N, self.calc_no_P, self.calc_p_N_edge, self.calc_n_P_edge,
                                            self.calc_N_B, self.calc_V_B, self.calc_mu_pP, self.calc_mu_nN, self.calc_mu_pN, self.calc_mu_nP,
                                            self.calc_D_nP, self.calc_D_pN, self.calc_D_nN, self.calc_D_pP,
                                            self.calc_tau_recP, self.calc_tau_recN, self.calc_tau_genP, self.calc_tau_genN,
                                            self.calc_L_pN, self.calc_L_nP, self.calc_J_nxdiffP, self.calc_J_pxdiffN])

        diode_quantities = {
            'N_a': None, #Acceptor concentration
            'N_d': None, #Donor concentration
            'V_bi': None, #Built-in voltage
            'V_bin': None, #Built-in voltage for n-type side
            'V_bip': None, #Built-in voltage for p-type side
            'W_dN': None, #Depletion width for n-type side
            'W_dP': None, #Depletion width for p-type side
            'V_PN': None, #PN junction voltage
            'po_N': None, #Hole concentration in n-type side at thermal equilibrium (minority carrier)
            'no_P': None, #Electron concentration in p-type side at thermal equilibrium (minority carrier)
            'p_N_edge': None, #Hole concentration at the edge of the depletion region in n-type side (minority carrier)
            'n_P_edge': None, #Electron concentration at the edge of the depletion region in p-type side (minority carrier)
            'V_B' : None, #Breakdown voltage
            'N_B' : None, #Breakdown concentration (more lightly doped concentration)
            'V_B' : None, #Breakdown voltage
            'mu_nP' : None, #Electron mobility in p-type side
            'mu_pP' : None, #Hole mobility in p-type side
            'mu_nN' : None, #Electron mobility in n-type side
            'mu_pN' : None, #Hole mobility in n-type side
            'D_nP' : None, #Diffusion coefficient for electrons in p-type side
            'D_pN' : None, #Diffusion coefficient for holes in n-type side
            'D_nN' : None, #Diffusion coefficient for electrons in n-type side
            'D_pP' : None, #Diffusion coefficient for holes in p-type side
            'tau_recP' : None, #Recombination time for p-type side
            'tau_recN' : None, #Recombination time for n-type side
            'tau_genP' : None, #Generation time for p-type side
            'tau_genN' : None, #Generation time for n-type side
            'L_pN' : None, #Diffusion length for holes on n-type side
            'L_nP' : None, #Diffusion length for electrons on p-type side
            'J_nxdiffP' : None, #Electronurrent density for p-type side diffusion current
            'J_pxdiffN' : None #Hole current density for n-type side diffusion current
        }

        self.known_quantities.update(diode_quantities)


    def calc_V_bi(self):
        if not self.should_calculate(needed=['T', 'N_a', 'N_d']):
            return
        if self.known_quantities['T'] != 300 * kelvin:
            return
        T = self.known_quantities['T']
        N_a = self.known_quantities['N_a']
        N_d = self.known_quantities['N_d']
        v_bi = k_BT_300K_DIVIDED_q * sp.log(N_a * N_d / n_i_300K**2)
        self.known_quantities['V_bi'] = v_bi

    def calc_V_bin(self):
        if not self.should_calculate(needed=['N_a', 'N_d', 'V_bi']):
            return
        V_bi = self.known_quantities['V_bi']
        N_a = self.known_quantities['N_a']
        N_d = self.known_quantities['N_d']
        v_bin = V_bi * N_d/(N_a + N_d)
        self.known_quantities['V_bin'] = v_bin

    def calc_V_bip(self):
        if not self.should_calculate(needed=['N_a', 'N_d', 'V_bi']):
            return
        V_bi = self.known_quantities['V_bi']
        N_a = self.known_quantities['N_a']
        N_d = self.known_quantities['N_d']
        v_bip = V_bi * N_a/(N_a + N_d)
        self.known_quantities['V_bip'] = v_bip

    def calc_W_dP(self):
        if not self.should_calculate(needed=['N_a', 'N_d', 'V_bi', 'V_PN']):
            return
        V_bi = remove_units(self.known_quantities['V_bi'])
        N_a = remove_units(self.known_quantities['N_a'])
        N_d = remove_units(self.known_quantities['N_d'])
        V_PN = remove_units(self.known_quantities['V_PN'])
        W_dp = sp.sqrt((2 * epsilon_Si_ / q_) * (N_d / (N_a * (N_a + N_d))) * (V_bi - V_PN)) * cm
        self.known_quantities['W_dP'] = W_dp

    def calc_W_dN(self):
        if not self.should_calculate(needed=['N_a', 'N_d', 'V_bi', 'V_PN']):
            return
        V_bi = remove_units(self.known_quantities['V_bi'])
        N_a = remove_units(self.known_quantities['N_a'])
        N_d = remove_units(self.known_quantities['N_d'])
        V_PN = remove_units(self.known_quantities['V_PN'])
        W_dn = sp.sqrt((2 * epsilon_Si_ / q_) * (N_a / (N_d * (N_a + N_d))) * (V_bi - V_PN)) * cm
        self.known_quantities['W_dN'] = W_dn

    def calc_po_N(self):
        if not self.should_calculate(needed=['N_d', 'T']):
            return
        if self.known_quantities['T'] != 300 * kelvin:
            return
        N_d = self.known_quantities['N_d']
        po_n = n_i_300K**2 / N_d
        self.known_quantities['po_N'] = po_n

    def calc_no_P(self):
        if not self.should_calculate(needed=['N_a', 'T']):
            return
        if self.known_quantities['T'] != 300 * kelvin:
            return
        N_a = self.known_quantities['N_a']
        no_p = n_i_300K**2 / N_a
        self.known_quantities['no_P'] = no_p

    def calc_p_N_edge(self):
        if not self.should_calculate(needed=['po_N', 'V_PN', 'T']):
            return
        po_N = self.known_quantities['po_N']
        p_n_edge = po_N * sp.exp(q * self.known_quantities['V_PN'] / (k_B * self.known_quantities['T']) / ELECTRON_VOLTS_PER_COLOUMB_TO_VOLTS)
        self.known_quantities['p_N_edge'] = p_n_edge

    def calc_n_P_edge(self):
        if not self.should_calculate(needed=['no_P', 'V_PN', 'T']):
            return
        no_P = self.known_quantities['no_P']
        n_p_edge = no_P * sp.exp(q * self.known_quantities['V_PN'] / (k_B * self.known_quantities['T']) / ELECTRON_VOLTS_PER_COLOUMB_TO_VOLTS)
        self.known_quantities['n_P_edge'] = n_p_edge

    def calc_N_B(self):
        if not self.should_calculate(needed=['N_a', 'N_d']):
            return
        N_a = self.known_quantities['N_a']
        N_d = self.known_quantities['N_d']
        N_B = min(N_a, N_d)
        self.known_quantities['N_B'] = N_B

    def calc_V_B(self):
        if not self.should_calculate(needed=['N_B']):
            return
        N_B = remove_units(self.known_quantities['N_B'])
        V_B = 60 * (E_gSi_300K_ / 1.1242) ** (3/2) * (1e16 / N_B) ** (3/4) * volt
        self.known_quantities['V_B'] = V_B

    def calc_mu_pP(self):
        if not self.should_calculate(needed=['N_a']):
            return
        N_a = remove_units(self.known_quantities['N_a'])
        mu_pP = (49.7 + 418.3 / (1+(N_a/(1.6e17))**(0.70))) * cm**2 / (volt * second)
        self.known_quantities['mu_pP'] = mu_pP

    def calc_mu_nP(self):
        if not self.should_calculate(needed=['N_a']):
            return
        N_a = remove_units(self.known_quantities['N_a'])
        mu_nP = (232 + 1180 / (1+(N_a/(8.0e16))**(0.90))) * cm**2 / (volt * second)
        self.known_quantities['mu_nP'] = mu_nP

    def calc_mu_nN(self):
        if not self.should_calculate(needed=['N_d']):
            return
        N_d = remove_units(self.known_quantities['N_d'])
        mu_nN = (92 + 1268 / (1+(N_d/(1.3e17))**(0.91))) * cm**2 / (volt * second)
        self.known_quantities['mu_nN'] = mu_nN

    def calc_mu_pN(self):
        if not self.should_calculate(needed=['N_d']):
            return
        N_d = remove_units(self.known_quantities['N_d'])
        mu_pN = (130 + 370 / (1+(N_d/(8.0e17))**(1.25))) * cm**2 / (volt * second)
        self.known_quantities['mu_pN'] = mu_pN

    def calc_D_nP(self):
        if not self.should_calculate(needed=['mu_nP', 'T']):
            return
        if self.known_quantities['T'] != 300 * kelvin:
            return
        mu_nP = self.known_quantities['mu_nP']
        D_nP = mu_nP * k_BT_300K_DIVIDED_q
        self.known_quantities['D_nP'] = D_nP

    def calc_D_pP(self):
        if not self.should_calculate(needed=['mu_pP', 'T']):
            return
        if self.known_quantities['T'] != 300 * kelvin:
            return
        mu_pP = self.known_quantities['mu_pP']
        D_pP = mu_pP * k_BT_300K_DIVIDED_q
        self.known_quantities['D_pP'] = D_pP

    def calc_D_nN(self):
        if not self.should_calculate(needed=['mu_nN', 'T']):
            return
        if self.known_quantities['T'] != 300 * kelvin:
            return
        mu_nN = self.known_quantities['mu_nN']
        D_nN = mu_nN * k_BT_300K_DIVIDED_q
        self.known_quantities['D_nN'] = D_nN

    def calc_D_pN(self):
        if not self.should_calculate(needed=['mu_pN', 'T']):
            return
        if self.known_quantities['T'] != 300 * kelvin:
            return
        mu_pN = self.known_quantities['mu_pN']
        D_pN = mu_pN * k_BT_300K_DIVIDED_q
        self.known_quantities['D_pN'] = D_pN

    def calc_tau_recP(self):
        if not self.should_calculate(needed=['N_a']):
            return
        N_a = remove_units(self.known_quantities['N_a'])
        tau_recP = 1 / (3.45e-12 * N_a + 9.50e-32 * (N_a**2)) * seconds
        self.known_quantities['tau_recP'] = tau_recP

    def calc_tau_recN(self):
        if not self.should_calculate(needed=['N_d']):
            return
        N_d = remove_units(self.known_quantities['N_d'])
        tau_recN = 1 / (7.80e-13 * N_d + 1.80e-31 * (N_d**2)) * seconds
        self.known_quantities['tau_recN'] = tau_recN

    def calc_tau_genP(self):
        if not self.should_calculate(needed=['tau_recP']):
            return
        tau_recP = self.known_quantities['tau_recP']
        tau_genP = 75 * tau_recP
        self.known_quantities['tau_genP'] = tau_genP

    def calc_tau_genN(self):
        if not self.should_calculate(needed=['tau_recN']):
            return
        tau_recN = self.known_quantities['tau_recN']
        tau_genN = 75 * tau_recN
        self.known_quantities['tau_genN'] = tau_genN

    def calc_L_pN(self):
        if not self.should_calculate(needed=['V_PN', 'tau_genN', 'tau_recN', 'D_pN']):
            return
        V_PN = self.known_quantities['V_PN']
        tau_genN = self.known_quantities['tau_genN']
        tau_recN = self.known_quantities['tau_recN']
        D_pN = self.known_quantities['D_pN']
        L_pN = None
        if V_PN > 0:
            L_pN = sp.sqrt(tau_recN * D_pN)
        else:
            L_pN = sp.sqrt(tau_genN * D_pN)
        self.known_quantities['L_pN'] = L_pN

    def calc_L_nP(self):
        if not self.should_calculate(needed=['V_PN', 'tau_genP', 'tau_recP', 'D_nP']):
            return
        V_PN = self.known_quantities['V_PN']
        tau_genP = self.known_quantities['tau_genP']
        tau_recP = self.known_quantities['tau_recP']
        D_nP = self.known_quantities['D_nP']
        L_nP = None
        if V_PN > 0:
            L_nP = sp.sqrt(tau_recP * D_nP)
        else:
            L_nP = sp.sqrt(tau_genP * D_nP)
        self.known_quantities['L_nP'] = L_nP

    def calc_J_pxdiffN(self):
        if not self.should_calculate(needed=['D_pN', 'L_pN', 'N_d', 'V_PN', 'T']):
            return
        if self.known_quantities['T'] != 300 * kelvin:
            return
        D_pN = self.known_quantities['D_pN']
        L_pN = self.known_quantities['L_pN']
        N_d = self.known_quantities['N_d']
        V_PN = self.known_quantities['V_PN']
        J_pxdiffN = (q * D_pN * n_i_300K**2) / (L_pN * N_d) * (sp.exp(V_PN / k_BT_300K_DIVIDED_q) - 1) * COULOUMBS_PER_SECOND_TO_AMPERES
        self.known_quantities['J_pxdiffN'] = J_pxdiffN

    def calc_J_nxdiffP(self):
        if not self.should_calculate(needed=['D_nP', 'L_nP', 'N_a', 'V_PN', 'T']):
            return
        if self.known_quantities['T'] != 300 * kelvin:
            return
        D_nP = self.known_quantities['D_nP']
        L_nP = self.known_quantities['L_nP']
        N_a = self.known_quantities['N_a']
        V_PN = self.known_quantities['V_PN']
        J_nxdiffP = (q * D_nP * n_i_300K**2) / (L_nP * N_a) * (sp.exp(V_PN / k_BT_300K_DIVIDED_q) - 1) * COULOUMBS_PER_SECOND_TO_AMPERES
        self.known_quantities['J_nxdiffP'] = J_nxdiffP