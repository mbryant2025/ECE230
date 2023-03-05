from ece230 import ECE230
from constants import *
from sympy.physics.units import *


class Silicon(ECE230):

    def __init__(self, T = 300 * kelvin):
        super().__init__(T)
        self.calculation_functions.extend([self.calc_Si_type, self.calc_n_oN, self.calc_p_oN, self.calc_n_oP, self.calc_p_oP])

        silicon_quantities = {
            'Si_type': None, #Type of silicon (n or p, or intrinsic)
            'N_a': None, #Acceptor concentration
            'N_d': None, #Donor concentration
            'n_oI': None, #Intrinsic electron concentration
            'p_oI': None, #Intrinsic hole concentration
            'n_oN': None, #Thermal equilibrium electron concentration in n-type silicon
            'p_oN': None, #Thermal equilibrium hole concentration in n-type silicon
            'n_oP': None, #Thermal equilibrium electron concentration in p-type silicon
            'p_oP': None, #Thermal equilibrium hole concentration in p-type silicon

        }

        self.known_quantities.update(silicon_quantities)

    def calc_Si_type(self):
        #if both Na and Nd are known, then the type of silicon is known as the greater
        #if only one is known, then the type of silicon is known as the known one
        #if neither is known, then the type of silicon is intrinsic
        if self.should_calculate(needed=['N_a', 'N_d']):
            if remove_units(self.get_known_quantity('N_a')) > remove_units(self.get_known_quantity('N_d')):
                self.known_quantities['Si_type'] = 'p'
            else:
                self.known_quantities['Si_type'] = 'n'
        elif self.should_calculate(needed=['N_a']):
            self.known_quantities['Si_type'] = 'p'
        elif self.should_calculate(needed=['N_d']):
            self.known_quantities['Si_type'] = 'n'
        else:
            self.known_quantities['Si_type'] = 'intrinsic'

    def calc_p_oP(self):
       #if N_a is known, then p_oP is known as N_a
        if self.should_calculate(needed=['N_a']):
            self.known_quantities['p_oP'] = self.known_quantities['N_a']
        #if n_oP is known, then p_oP is known as n_i**2 / n_oP
        elif self.should_calculate(needed=['n_oP']):
            self.known_quantities['p_oP'] = n_i_300K**2 / self.known_quantities['n_oP']

    def calc_n_oP(self):
        #if p_oP is known, then n_oP is known as ni**2 / p_oP
        if self.should_calculate(needed=['p_oP']):
            self.known_quantities['n_oP'] = n_i_300K**2 / self.known_quantities['p_oP']

    
    def calc_n_oN(self):
        #if N_d is known, then n_oN is known as N_d
        if self.should_calculate(needed=['N_d']):
            self.known_quantities['n_oN'] = self.known_quantities['N_d']
        #if p_oN is known, then n_oN is known as n_i**2 / p_oN
        elif self.should_calculate(needed=['p_oN']):
            self.known_quantities['n_oN'] = n_i_300K**2 / self.known_quantities['p_oN']

    def calc_p_oN(self):
        #if n_oN is known, then p_oN is known as n_i**2 / n_oN
        if self.should_calculate(needed=['n_oN']):
            self.known_quantities['p_oN'] = n_i_300K**2 / self.known_quantities['n_oN']