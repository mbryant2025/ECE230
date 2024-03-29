from silicon import Silicon
from constants import *
import sympy as sp
from sympy.physics.units import *


class Diode(Silicon):

    def __init__(self, T = 300 * kelvin):
        super().__init__(T)
        self.calculation_functions.extend([self.calc_bias_mode, self.calc_V_bi, self.calc_V_bin, self.calc_V_bip, self.calc_W_d, self.calc_W_dN, self.calc_W_dP,
                                           self.calc_rho_p, self.calc_rho_n, self.calc_epsilon_x_p, self.calc_epsilon_x_n,
                                            self.calc_po_N, self.calc_no_P, self.calc_p_N_edge, self.calc_n_P_edge,
                                            self.calc_N_B, self.calc_V_B, self.calc_mu_pP, self.calc_mu_nN, self.calc_mu_pN, self.calc_mu_nP,
                                            self.calc_D_nP, self.calc_D_pN, self.calc_D_nN, self.calc_D_pP,
                                            self.calc_tau_recP, self.calc_tau_recN, self.calc_tau_genP, self.calc_tau_genN,
                                            self.calc_L_pN, self.calc_L_nP, self.calc_pN, self.calc_nP, self.calc_J_nxdiffP, self.calc_J_pxdiffN,
                                            self.calc_C_pndep_per_area, self.calc_C_pndep, self.calc_epsilon_x_max,
                                            self.calc_tau_recSCR, self.calc_J_Dscr, self.calc_J_Sdiff, self.calc_J_Sscr, self.calc_J_D, self.calc_I_D,
                                            self.calc_tau_genSCR, self.calc_J_Dscg, self.calc_Q_pndep, self.calc_Q_pndiffP, self.calc_Q_pndiffN,
                                            self.calc_C_pndiffN, self.calc_C_pndiffP, self.calc_C_pndiff, self.calc_C_pn, self.calc_L_p_rec_N, self.calc_L_p_gen_N, self.calc_N_d,
                                           self.calc_N_a, self.calc_L_n_rec_P, self.calc_L_n_gen_P, self.calc_J_Ddiff, self.calc_I_Ddiff
                                           ])

        diode_quantities = {
            'Bias_Mode': None, #Bias mode (forward or reverse)
            'N_a': None, #Acceptor concentration
            'N_d': None, #Donor concentration
            'V_bi': None, #Built-in voltage
            'V_bin': None, #Built-in voltage for n-type side
            'V_bip': None, #Built-in voltage for p-type side
            'W_dN': None, #Depletion width for n-type side
            'W_dP': None, #Depletion width for p-type side
            'W_d': None, #Depletion width
            'W_N': None, #Width of n-type side NOTE user input
            'W_P': None, #Width of p-type side NOTE user input
            'V_PN': 0 * volt, #PN junction voltage #NOTE user input
            'rho_p' : None, #net-charge-concentration distribution in p side depletion region
            'rho_n' : None, #net-charge-concentration distribution in n side depletion region
            'epsilon_x_max' : None, #Maximum e-field in depletion region
            'epsilon_x_p' : None, #e-field in p-type side depletion region
            'epsilon_x_n' : None, #e-field in n-type side depletion region
            'po_N': None, #Hole concentration in n-type side at thermal equilibrium (minority carrier)
            'no_P': None, #Electron concentration in p-type side at thermal equilibrium (minority carrier)
            'pN' : None, #The hole-concentration distribution in a long-base n-side quasi-neutral region
            'nP' : None, #The electron concentration distribution in a long-base p-side quasi-neutral region
            'x' : None, #The x coordinate NOTE: user input
            'p_N_edge': None, #Hole concentration at the edge of the depletion region in n-type side (minority carrier)
            'n_P_edge': None, #Electron concentration at the edge of the depletion region in p-type side (minority carrier)
            'V_B' : None, #Breakdown voltage
            'N_B' : None, #Breakdown concentration (more lightly doped concentration)
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
            'L_p_rec_N' : None, #Recombination length for holes on n-type side
            'L_p_gen_N' : None, #Generation length for electrons on p-type side
            'L_n_rec_P' : None, #Recombination length for electrons on p-type side
            'J_nxdiffP' : None, #Electronurrent density for p-type side diffusion current
            'J_pxdiffN' : None, #Hole current density for n-type side diffusion current
            'J_Ddiff' : None, #Diode diffusion current density
            'I_Ddiff' : None, #Diode diffusion current
            'C_pn*dep_per_area' : None, #Capacitance of depletion region per unit area
            'C_pn*dep' : None, #Capacitance of depletion region
            'A' : None, #Area NOTE: user input
            'tau_recSCR' : None, #Space-charge-region recombination lifetime
            'tau_genSCR' : None, #Space-charge-region generation lifetime
            'J_Dscr' : None, #Diode space-charge-recombination current
            'J_Sdiff' : None, #Diffusion saturation current density
            'J_Sscr' : None, #Space-charge-recombination saturation current density
            'J_D' : None, #Diode current density
            'I_D' : None, #Diode current
            'J_Dscg' : None, #Diode space-charge-generation current density (reverse bias)
            'Q_pndep' : None, #Depletion charge density
            'Q_pndiffP' : None, #p-side quais-neutral region electron diffusion charge density
            'Q_pndiffN' : None, #n-side quais-neutral region hole diffusion charge density
            'C_pn*diffN' : None, #Small-signal diffusion capacitance per unit area onthe n-side of the junction
            'C_pn*diffP' : None, #Small-signal diffusion capacitance per unit area onthe p-side of the junction
            'C_pn*diff' : None, #Small-signal diffusion capacitance per unit area
            'C_pn' : None #Small-signal capacitance per unit area
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

    def calc_bias_mode(self):
        if not self.should_calculate(needed=['V_PN']):
            return
        V_PN = self.known_quantities['V_PN']
        if V_PN >= 0 * volt:
            self.known_quantities['Bias_Mode'] = 'Forward'
        else:
            self.known_quantities['Bias_Mode'] = 'Reverse'

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

    def calc_W_d(self):
        if self.should_calculate(needed=['N_a', 'N_d', 'V_bi', 'V_PN']):
            V_PN = self.known_quantities['V_PN']
            V_bi = self.known_quantities['V_bi']
            N_a = self.known_quantities['N_a']
            N_d = self.known_quantities['N_d']
            W_d = remove_units(sp.sqrt((2 * epsilon_Si / q) * (N_a + N_d) / (N_a * N_d) * (V_bi-V_PN))) * cm
            self.known_quantities['W_d'] = W_d

    def calc_W_dP(self):
        if not self.should_calculate(needed=['W_d', 'N_a', 'N_d']):
            return
        W_d = self.known_quantities['W_d']
        N_a = self.known_quantities['N_a']
        N_d = self.known_quantities['N_d']
        W_dP = W_d * N_d/(N_a + N_d)
        self.known_quantities['W_dP'] = W_dP

    def calc_W_dN(self):
        if not self.should_calculate(needed=['W_d', 'N_a', 'N_d']):
            return
        W_d = self.known_quantities['W_d']
        N_a = self.known_quantities['N_a']
        N_d = self.known_quantities['N_d']
        W_dN = W_d * N_a/(N_a + N_d)
        self.known_quantities['W_dN'] = W_dN

    def calc_rho_p(self):
        if not self.should_calculate(needed=['N_a']):
            return
        N_a = self.known_quantities['N_a']
        rho_p = -q * N_a
        self.known_quantities['rho_p'] = rho_p

    def calc_rho_n(self):
        if not self.should_calculate(needed=['N_d']):
            return
        N_d = self.known_quantities['N_d']
        rho_n = q * N_d
        self.known_quantities['rho_n'] = rho_n

    def calc_epsilon_x_p(self):
        if not self.should_calculate(needed=['N_a', 'W_dP', 'x']):
            return
        N_a = self.known_quantities['N_a']
        W_dP = self.known_quantities['W_dP']
        x = self.known_quantities['x']
        epsilon_x_p = -q * N_a * epsilon_Si * (x + W_dP)
        self.known_quantities['epsilon_x_p'] = epsilon_x_p

    def calc_epsilon_x_n(self):
        if not self.should_calculate(needed=['N_d', 'W_dN', 'x']):
            return
        N_d = self.known_quantities['N_d']
        W_dN = self.known_quantities['W_dN']
        x = self.known_quantities['x']
        epsilon_x_n = -q * N_d * epsilon_Si * (W_dN - x)
        self.known_quantities['epsilon_x_n'] = epsilon_x_n

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

    def calc_N_d(self):
        if self.should_calculate(needed=['mu_nN']):
            mu_nN = remove_units(self.known_quantities['mu_nN'])
            N_d = 1.3e17 * ((1268/(mu_nN-92))-1) ** (1/0.91) * cm**-3
            self.known_quantities['N_d'] = N_d
        elif self.should_calculate(needed=['mu_pN']):
            mu_pN = remove_units(self.known_quantities['mu_pN'])
            N_d = 8.0e17 * ((370/(mu_pN-130))-1) ** (1/1.25) * cm**-3
            self.known_quantities['N_d'] = N_d

    def calc_N_a(self):
        if self.should_calculate(needed=['mu_nP']):
            mu_pP = remove_units(self.known_quantities['mu_nP'])
            N_a = 8.0e16 * ((1180/(mu_pP-232))-1) ** (1/0.90) * cm**-3
            self.known_quantities['N_a'] = N_a
        elif self.should_calculate(needed=['mu_pP']):
            mu_pP = remove_units(self.known_quantities['mu_pP'])
            N_a = 1.6e17 * ((418.3/(mu_pP-49.7))-1) ** (1/0.70) * cm**-3
            self.known_quantities['N_a'] = N_a

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

    def calc_L_p_rec_N(self):
        if not self.should_calculate(needed=['tau_recN', 'D_pN']):
            return
        tau_recN = self.known_quantities['tau_recN']
        D_pN = self.known_quantities['D_pN']
        L_p_rec_N = sp.sqrt(tau_recN * D_pN)
        self.known_quantities['L_p_rec_N'] = L_p_rec_N

    def calc_L_p_gen_N(self):
        if not self.should_calculate(needed=['tau_genN', 'D_pN']):
            return
        tau_genN = self.known_quantities['tau_genN']
        D_pN = self.known_quantities['D_pN']
        L_p_gen_N = sp.sqrt(tau_genN * D_pN)
        self.known_quantities['L_p_gen_N'] = L_p_gen_N

    def calc_L_n_rec_P(self):
        if not self.should_calculate(needed=['tau_recP', 'D_nP']):
            return
        tau_recP = self.known_quantities['tau_recP']
        D_nP = self.known_quantities['D_nP']
        L_n_rec_P = sp.sqrt(tau_recP * D_nP)
        self.known_quantities['L_n_rec_P'] = L_n_rec_P

    def calc_L_n_gen_P(self):
        if not self.should_calculate(needed=['tau_genP', 'D_nP']):
            return
        tau_genP = self.known_quantities['tau_genP']
        D_nP = self.known_quantities['D_nP']
        L_n_gen_P = sp.sqrt(tau_genP * D_nP)
        self.known_quantities['L_n_gen_P'] = L_n_gen_P

    def calc_pN(self):
        if not self.should_calculate(needed=['V_PN', 'x', 'W_dN', 'L_pN', 'po_N']):
            return
        x = self.known_quantities['x']
        V_PN = self.known_quantities['V_PN']
        W_dN = self.known_quantities['W_dN']
        L_pN = self.known_quantities['L_pN']
        po_N = self.known_quantities['po_N']
        pN = po_N * sp.exp(-(x - W_dN) / L_pN) * (sp.exp(V_PN / k_BT_300K_DIVIDED_q) - 1) + po_N
        self.known_quantities['pN'] = pN
    
    def _calc_pN_for_x(self, x):
        if not self.should_calculate(needed=['V_PN', 'W_dN', 'L_pN', 'po_N']):
            return 0
        V_PN = self.known_quantities['V_PN']
        W_dN = self.known_quantities['W_dN']
        L_pN = self.known_quantities['L_pN']
        po_N = self.known_quantities['po_N']
        pN = po_N * sp.exp(-(x - W_dN) / L_pN) * (sp.exp(V_PN / k_BT_300K_DIVIDED_q) - 1) + po_N
        return remove_units(pN)

    def calc_nP(self):
        if not self.should_calculate(needed=['V_PN', 'x', 'W_dP', 'L_nP', 'no_P']):
            return
        x = self.known_quantities['x']
        V_PN = self.known_quantities['V_PN']
        W_dP = self.known_quantities['W_dP']
        L_nP = self.known_quantities['L_nP']
        no_P = self.known_quantities['no_P']
        nP = no_P * sp.exp((x + W_dP) / L_nP) * (sp.exp(V_PN / k_BT_300K_DIVIDED_q) - 1) + no_P
        self.known_quantities['nP'] = nP

    def _calc_nP_for_x(self, x):
        if not self.should_calculate(needed=['V_PN', 'W_dP', 'L_nP', 'no_P']):
            return 0
        V_PN = self.known_quantities['V_PN']
        W_dP = self.known_quantities['W_dP']
        L_nP = self.known_quantities['L_nP']
        no_P = self.known_quantities['no_P']
        nP = no_P * sp.exp((x + W_dP) / L_nP) * (sp.exp(V_PN / k_BT_300K_DIVIDED_q) - 1) + no_P
        return remove_units(nP)

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
        W_dP = self.known_quantities['W_dP']
        J_nxdiffP = (q * D_nP * n_i_300K**2) / (L_nP * N_a) * (sp.exp(V_PN / k_BT_300K_DIVIDED_q) - 1) * COULOUMBS_PER_SECOND_TO_AMPERES
        self.known_quantities['J_nxdiffP'] = J_nxdiffP

    def calc_I_Ddiff(self):
        if not self.should_calculate(needed=['J_pxdiffN', 'J_nxdiffP', 'A']):
            return
        J_pxdiffN = self.known_quantities['J_pxdiffN']
        J_nxdiffP = self.known_quantities['J_nxdiffP']
        A = self.known_quantities['A']
        I_Ddiff = (J_pxdiffN + J_nxdiffP) * A
        self.known_quantities['I_Ddiff'] = I_Ddiff

    def calc_J_Ddiff(self):
        if not self.should_calculate(needed=['J_pxdiffN', 'J_nxdiffP']):
            return
        J_pxdiffN = self.known_quantities['J_pxdiffN']
        J_nxdiffP = self.known_quantities['J_nxdiffP']
        J_Ddiff = J_pxdiffN + J_nxdiffP
        self.known_quantities['J_Ddiff'] = J_Ddiff

    def calc_C_pndep_per_area(self):
        if not self.should_calculate(needed=['N_d', 'N_a', 'V_PN', 'V_bi']):
            return
        N_a = remove_units(self.known_quantities['N_a'])
        N_d = remove_units(self.known_quantities['N_d'])
        V_PN = remove_units(self.known_quantities['V_PN'])
        V_bi = remove_units(self.known_quantities['V_bi'])
        C_pndep_per_area = sp.sqrt((q_ * epsilon_Si_ * N_a * N_d) / (2 * (N_d + N_a) * (V_bi - V_PN))) * farad / cm**2
        self.known_quantities['C_pn*dep_per_area'] = C_pndep_per_area

    def calc_C_pndep(self):
        if not self.should_calculate(needed=['C_pn*dep_per_area', 'A']):
            return
        A = self.known_quantities['A']
        C_pndep_per_area = self.known_quantities['C_pn*dep_per_area']
        C_pndep = C_pndep_per_area * A
        self.known_quantities['C_pn*dep'] = C_pndep

    def calc_epsilon_x_max(self):
        if not self.should_calculate(needed=['N_d', 'W_dN']):
            return
        N_d = self.known_quantities['N_d']
        W_dN = self.known_quantities['W_dN']
        epsilon_x_max = -q * N_d * W_dN / epsilon_Si
        self.known_quantities['epsilon_x_max'] = epsilon_x_max

    def calc_tau_recSCR(self):
        if not self.should_calculate(needed=['tau_recN', 'tau_recP']):
            return
        tau_recN = self.known_quantities['tau_recN']
        tau_recP = self.known_quantities['tau_recP']
        tau_recSCR = (tau_recN + tau_recP) / 2
        self.known_quantities['tau_recSCR'] = tau_recSCR

    def calc_J_Dscr(self):
        if not self.should_calculate(needed=['W_d', 'tau_recSCR', 'V_PN']):
            return
        W_d = self.known_quantities['W_d']
        tau_recSCR = self.known_quantities['tau_recSCR']
        V_PN = self.known_quantities['V_PN']
        J_Dscr = (q * n_i_300K * W_d) / (2 * tau_recSCR) * (sp.exp(V_PN / (2 * k_BT_300K_DIVIDED_q)) - 1) * COULOUMBS_PER_SECOND_TO_AMPERES
        self.known_quantities['J_Dscr'] = J_Dscr

    def calc_J_Sdiff(self):
        if not self.should_calculate(needed=['D_nP', 'L_nP', 'N_a', 'D_pN', 'L_pN', 'N_d']):
            return
        D_nP = self.known_quantities['D_nP']
        L_nP = self.known_quantities['L_nP']
        N_a = self.known_quantities['N_a']
        D_pN = self.known_quantities['D_pN']
        L_pN = self.known_quantities['L_pN']
        N_d = self.known_quantities['N_d']
        J_Sdiff = (q * D_nP * n_i_300K**2) / (L_nP * N_a) + (q * D_pN * n_i_300K**2) / (L_pN * N_d)
        self.known_quantities['J_Sdiff'] = J_Sdiff

    def calc_J_Sscr(self):
        if not self.should_calculate(needed=['W_d', 'tau_recSCR']):
            return
        W_d = self.known_quantities['W_d']
        tau_recSCR = self.known_quantities['tau_recSCR']
        J_Sscr = (q * n_i_300K * W_d) / (2 * tau_recSCR)
        self.known_quantities['J_Sscr'] = J_Sscr

    def calc_J_D(self):
        if not self.should_calculate(needed=['J_Sdiff', 'J_Sscr', 'V_PN']):
            return
        J_Sdiff = self.known_quantities['J_Sdiff']
        J_Sscr = self.known_quantities['J_Sscr']
        V_PN = self.known_quantities['V_PN']
        J_D = J_Sdiff * (sp.exp(V_PN / (k_BT_300K_DIVIDED_q)) - 1) + J_Sscr * (sp.exp(V_PN / (2 * k_BT_300K_DIVIDED_q)) - 1)
        self.known_quantities['J_D'] = J_D

    def calc_I_D(self):
        if not self.should_calculate(needed=['J_D', 'A']):
            return
        J_D = self.known_quantities['J_D']
        A = self.known_quantities['A']
        I_D = J_D * A
        self.known_quantities['I_D'] = I_D * COULOUMBS_PER_SECOND_TO_AMPERES

    def calc_tau_genSCR(self):
        #factor of 75 as before
        if not self.should_calculate(needed=['tau_recSCR']):
            return
        tau_recSCR = self.known_quantities['tau_recSCR']
        tau_genSCR = tau_recSCR * 75
        self.known_quantities['tau_genSCR'] = tau_genSCR


    def calc_J_Dscg(self):
        if not self.should_calculate(needed=['W_d', 'tau_genSCR', 'V_PN']):
            return
        W_d = self.known_quantities['W_d']
        tau_recSCR = self.known_quantities['tau_genSCR']
        V_PN = self.known_quantities['V_PN']
        J_Dscg = (q * n_i_300K * W_d) / (2 * tau_recSCR) * (sp.exp(V_PN / (2 * k_BT_300K_DIVIDED_q)) - 1) * COULOUMBS_PER_SECOND_TO_AMPERES
        self.known_quantities['J_Dscg'] = J_Dscg

    def calc_Q_pndep(self):
        # if self.should_calculate(needed=['N_d', 'W_dN']):
        #     N_d = self.known_quantities['N_d']
        #     W_dN = self.known_quantities['W_dN']
        #     Q_pndep = q * N_d * W_dN
        #     self.known_quantities['Q_pndep'] = Q_pndep
        #     return
        # if self.should_calculate(needed=['N_a', 'W_dP']):
        #     N_a = self.known_quantities['N_a']
        #     W_dP = self.known_quantities['W_dP']
        #     Q_pndep = q * N_a * W_dP
        #     self.known_quantities['Q_pndep'] = Q_pndep
        #     return

        if not self.should_calculate(needed=['N_a', 'N_d', 'V_bi', 'V_PN']):
            return
        N_a = remove_units(self.known_quantities['N_a'])
        N_d = remove_units(self.known_quantities['N_d'])
        V_bi = remove_units(self.known_quantities['V_bi'])
        V_PN = remove_units(self.known_quantities['V_PN'])
        Q_pndep = sp.sqrt(2 * q_ * epsilon_Si_ * (N_a * N_d) / (N_a + N_d) * (V_bi - V_PN)) * coulombs / angstroms**2
        self.known_quantities['Q_pndep'] = Q_pndep
        
    def calc_Q_pndiffP(self):
        # if not self.should_calculate(needed=['no_P', 'W_P', 'W_dP']):
        #     return
        # no_p = remove_units(self.known_quantities['no_P'])
        # W_P = remove_units(self.known_quantities['W_P'])
        # W_dP = remove_units(self.known_quantities['W_dP'])
        # x = sp.symbols('x')
        # integrand = self._calc_nP_for_x(x) - no_p
        # integral_limits = (-W_P, -W_dP)
        # Q_pndiffP = -q * sp.integrate(integrand, (x, *integral_limits)) / angstroms**2
        # self.known_quantities['Q_pndiffP'] = Q_pndiffP
        #need L_nP, N_a, V_PN

        #The electron diffusion charge density in a long-base p-side quasi-neutral region is given by
        if not self.should_calculate(needed=['L_nP', 'N_a', 'V_PN']):
            return
        L_nP = self.known_quantities['L_nP']
        N_a = self.known_quantities['N_a']
        V_PN = self.known_quantities['V_PN']
        Q_pndiffP = -q * n_i_300K**2 * L_nP * (sp.exp(V_PN / k_BT_300K_DIVIDED_q) - 1) / N_a
        self.known_quantities['Q_pndiffP'] = Q_pndiffP

    def calc_Q_pndiffN(self):
        # if not self.should_calculate(needed=['po_N', 'W_N', 'W_dN']):
        #     return
        # po_n = remove_units(self.known_quantities['po_N'])
        # W_N = remove_units(self.known_quantities['W_N'])
        # W_dN = remove_units(self.known_quantities['W_dN'])
        # x = sp.symbols('x')
        # integrand = self._calc_pN_for_x(x) - po_n
        # integral_limits = (W_dN, W_N)
        # Q_pndiffN = q * sp.integrate(integrand, (x, *integral_limits)) / angstroms**2
        # self.known_quantities['Q_pndiffN'] = Q_pndiffN
        #need L_pN, N_d, V_PN

        #The hole diffusion charge density in a long-base n-side quasi-neutral region is given by
        if not self.should_calculate(needed=['L_pN', 'N_d', 'V_PN']):
            return
        L_pN = self.known_quantities['L_pN']
        N_d = self.known_quantities['N_d']
        V_PN = self.known_quantities['V_PN']
        Q_pndiffN = q * n_i_300K**2 * L_pN * (sp.exp(V_PN / k_BT_300K_DIVIDED_q) - 1) / N_d
        self.known_quantities['Q_pndiffN'] = Q_pndiffN

    def calc_C_pndiffN(self):
        if not self.should_calculate(needed=['L_p_rec_N', 'N_d', 'V_PN']):
            return
        L_pN = self.known_quantities['L_p_rec_N']
        N_d = self.known_quantities['N_d']
        V_PN = self.known_quantities['V_PN']
        C_pndiffN = remove_units(q**2 * n_i_300K**2 * L_pN * (sp.exp(V_PN / k_BT_300K_DIVIDED_q)) / (N_d * k_BT_300K)) * F / cm**2# * COULOUMBS_SQUARE_PER_EV_CENTIMETER_SQUARE_TO_FARADS_PER_ANGSTROM_SQUARE
        self.known_quantities['C_pn*diffN'] = C_pndiffN

    def calc_C_pndiffP(self):
        if not self.should_calculate(needed=['L_n_rec_P', 'N_a', 'V_PN']):
            return
        L_nP = self.known_quantities['L_n_rec_P']
        N_a = self.known_quantities['N_a']
        V_PN = self.known_quantities['V_PN']
        C_pndiffP = remove_units(q**2 * n_i_300K**2 * L_nP * (sp.exp(V_PN / k_BT_300K_DIVIDED_q)) / (N_a * k_BT_300K)) * F / cm**2 # * COULOUMBS_SQUARE_PER_EV_CENTIMETER_SQUARE_TO_FARADS_PER_ANGSTROM_SQUARE
        self.known_quantities['C_pn*diffP'] = C_pndiffP

    def calc_C_pndiff(self):
        if not self.should_calculate(needed=['C_pn*diffP', 'C_pn*diffN']):
            return
        C_pndiffP = self.known_quantities['C_pn*diffP']
        C_pndiffN = self.known_quantities['C_pn*diffN']
        C_pndiff = C_pndiffP + C_pndiffN
        self.known_quantities['C_pn*diff'] = C_pndiff

    def calc_C_pn(self):
        if not self.should_calculate(needed=['C_pn*diff', 'C_pn*dep_per_area']):
            return
        C_pndiff = self.known_quantities['C_pn*diff']
        C_pndep = self.known_quantities['C_pn*dep_per_area']
        C_pn = C_pndiff + C_pndep
        self.known_quantities['C_pn'] = C_pn