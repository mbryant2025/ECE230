from sympy.physics.units import *
from constants import *

class ECE230:

    def __init__(self, T = 300 * kelvin):
        self.calculation_functions = []
        self.known_quantities = {}
        self.known_quantities['T'] = T
    
    def check_if_known(self, quantities):
        if quantities is None:
            return False
        for quantity in quantities:
            if self.known_quantities[quantity] is None:
                return False
        return True
    
    def calculate_values(self):
        #repeats multiple times to ensure all values are calculated (dependencies)
        for _ in range(3):
            for calculation_function in self.calculation_functions:
                calculation_function()

    def get_known_quantities(self):
        return self.known_quantities
    
    def get_known_quantity(self, quantity):
        return self.known_quantities[quantity]
    
    def add_known_quantity(self, quantity, value):
        self.known_quantities[quantity] = value
        self.calculate_values()

    def print_known_quantities(self):
        self.calculate_values()
        for quantity in self.known_quantities:
            print(quantity, self.known_quantities[quantity])

    def get_possible_quantities(self):
        return self.known_quantities.keys()
    
    def should_calculate(self, needed=None):
        return self.check_if_known(needed)
