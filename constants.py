from sympy import pi
from sympy.physics import units
from sympy.physics.units import Quantity
import math

"""
This module contains constants and conversion factors

Constants with a trailing underscore are unitless.
Constants without a trailing underscore have corresponding units.

"""

def remove_units(expr):
    """Remove units from an expression.

    Parameters
    ----------
    expr : sympy expression

    Returns
    -------
    expr : sympy expression
        The expression with units removed.

    """
    if isinstance(expr, (int, float)):
        return expr
    if not expr.has(Quantity):
        return expr
    units = expr.subs({x: 1 for x in expr.args if not x.has(Quantity)})
    return expr / units

angstroms = units.cm / 1e8

a_Si = 5.43086 * angstroms
"""The lattice parameter of silicon in angstroms."""

a_Si_ = 5.43086
"""The lattice parameter of silicon in angstroms."""

N_Si = 5e22 * units.cm**-3
"""The density of silicon atoms in cm^-3."""

N_Si_ = 5e22
"""The density of silicon atoms in cm^-3."""

m_0 = 9.109389e-28 * units.g
"""The mass of an electron in g."""

m_0_ = 9.109389e-28
"""The mass of an electron in g."""

m_star_ndos_300K = 1.090 * m_0
"""The electron density of states effective mass at 300K in g."""

m_star_ndos_300K_ = 1.090 * m_0_
"""The electron density of states effective mass at 300K in g."""

m_star_ndos_400K = 1.110 * m_0
"""The electron density of states effective mass at 400K in g."""

m_star_ndos_400K_ = 1.110 * m_0_
"""The electron density of states effective mass at 400K in g."""

m_star_pdos_300K = 1.150 * m_0
"""The hole density of states effective mass at 300K in g."""

m_star_pdos_300K_ = 1.150 * m_0_
"""The hole density of states effective mass at 300K in g."""

m_star_pdos_400K = 1.230 * m_0
"""The hole density of states effective mass at 400K in g."""

m_star_pdos_400K_ = 1.230 * m_0_
"""The hole density of states effective mass at 400K in g."""

m_star_n = 0.2588 * m_0
"""The electron conductivity effective mass in g."""

m_star_n_ = 0.2588 * m_0_
"""The electron conductivity effective mass in g."""

m_star_p = 0.38158 * m_0
"""The hole effective mass in g."""

m_star_p_ = 0.38158 * m_0_
"""The hole effective mass in g."""

h = 4.135669e-15 * units.eV * units.s
"""The Planck constant in eV*s."""

h_ = 4.135669e-15
"""The Planck constant in eV*s."""

hbar = h / (2 * pi)
"""The reduced Planck constant in eV*s."""

hbar_ = h_ / (2 * math.pi)
"""The reduced Planck constant in eV*s."""

k_B = 8.617385e-5 * units.eV / units.K
"""The Boltzmann constant in eV/K."""

k_B_ = 8.617385e-5
"""The Boltzmann constant in eV/K."""

k_BT_300K = 4.1419479804e-21 * units.eV
"""The thermal energy at 300K in eV."""

k_BT_300K_ = 4.1419479804e-21
"""The thermal energy at 300K in eV."""

n_i_300K = 1.07e10 * units.cm**-3
"""The intrinsic carrier concentration at 300K in cm^-3."""

n_i_300K_ = 1.07e10
"""The intrinsic carrier concentration at 300K in cm^-3."""

N_c_300K = 2.86e19 * units.cm**-3
"""The conduction band effective density of states at 300K in cm^-3."""

N_c_300K_ = 2.86e19
"""The conduction band effective density of states at 300K in cm^-3."""

N_v_300K = 3.10e19 * units.cm**-3
"""The valence band effective density of states at 300K in cm^-3."""

N_v_300K_ = 3.10e19
"""The valence band effective density of states at 300K in cm^-3."""

ANGSTROMS_TO_METERS = 1e-10 * units.m / angstroms
"""The conversion factor from angstroms to meters."""

ANGSTROMS_TO_METERS_ = 1e-10
"""The conversion factor from angstroms to meters."""

ANGSTROMS_TO_CENTIMETERS = 1e-8 * units.cm / angstroms
"""The conversion factor from angstroms to centimeters."""

ANGSTROMS_TO_CENTIMETERS_ = 1e-8
"""The conversion factor from angstroms to centimeters."""

CENTIMETERS_TO_ANGSTROMS = 1e8 * angstroms / units.cm
"""The conversion factor from centimeters to angstroms."""

CENTIMETERS_TO_ANGSTROMS_ = 1e8
"""The conversion factor from centimeters to angstroms."""

CENTIMETERS_TO_METERS = 1e-2 * units.m / units.cm
"""The conversion factor from centimeters to meters."""

CENTIMETERS_TO_METERS_ = 1e-2
"""The conversion factor from centimeters to meters."""

METERS_TO_ANGSTROMS = 1e10 * angstroms / units.m
"""The conversion factor from meters to angstroms."""

METERS_TO_ANGSTROMS_ = 1e10
"""The conversion factor from meters to angstroms."""

q = 1.602177e-19 * units.C
"""The elementary charge in C."""

q_ = 1.602177e-19
"""The elementary charge in C."""

E_gSi_300K = 1.1242 * units.eV
"""The band gap of silicon at 300K in eV."""

E_gSi_300K_ = 1.1242
"""The band gap of silicon at 300K in eV."""

E_gSi_0k = 1.17 * units.eV
"""The band gap of silicon at 0K in eV."""

E_gSi_0k_ = 1.17
"""The band gap of silicon at 0K in eV."""

ELECTRON_VOLTS_PER_COLOUMB_TO_VOLTS = 1.602176487e-19 * units.volt * units.C / units.eV
"""The conversion factor from eV/C to V."""

ELECTRON_VOLTS_PER_COLOUMB_TO_VOLTS_ = 1.602176487e-19
"""The conversion factor from eV/C to V."""

epsilon0 = 8.854187e-14 * units.F / units.cm
"""The permittivity of free space in F/cm."""

epsilon0_ = 8.854187e-14
"""The permittivity of free space in F/cm."""

epsilon_Si = 11.7 * epsilon0
"""The permittivity of silicon in F/cm."""

epsilon_Si_ = 11.7 * epsilon0_
"""The permittivity of silicon in F/cm."""

epsilon_ox = 3.90 * epsilon0
"""The permittivity of oxide in F/cm."""

epsilon_ox_ = 3.90 * epsilon0_
"""The permittivity of oxide in F/cm."""

k_BT_300K_DIVIDED_q = 0.025852 * units.volt
"""The thermal energy at 300K in V."""

k_BT_300K_DIVIDED_q_ = 0.025852
"""The thermal energy at 300K in V."""

k_B_DIVIDED_q = k_B / q * ELECTRON_VOLTS_PER_COLOUMB_TO_VOLTS
"""The Boltzmann divided by the elementary charge in V/K."""

k_B_DIVIDED_q_ = k_B_ / q_ * ELECTRON_VOLTS_PER_COLOUMB_TO_VOLTS_
"""The Boltzmann divided by the elementary charge in V/K."""

root_mean_square_velocity_silicon_300K = 1.0e7 * units.cm / units.s
"""The root mean square thermal velocity of silicon in cm/s."""

root_mean_square_velocity_silicon_300K_ = 1.0e7
"""The root mean square thermal velocity of silicon in cm/s."""

v_n_sat = 1.066e7 * units.cm / units.s
"""The electron saturation velocity in cm/s."""

v_n_sat_ = 1.066e7
"""The electron saturation velocity in cm/s."""

v_p_sat = 8.290e6 * units.cm / units.s
"""The hole saturation velocity in cm/s."""

v_p_sat_ = 8.290e6
"""The hole saturation velocity in cm/s."""

COULOUMBS_PER_SECOND_TO_AMPERES = units.A / units.C * units.s
"""The conversion factor from C/s to A."""

COULOUMBS_SQUARE_PER_EV_CENTIMETER_SQUARE_TO_FARADS_PER_ANGSTROM_SQUARE = units.F / units.cm**2 * units.eV / units.C**2 * angstroms**2 * 6.24150974e18
"""The conversion factor from C^2/(eV cm^2) to F/angstrom^2."""

COULOUMB_PER_FARAD_TO_VOLT = units.volt / units.C * units.F
"""The conversion factor from C/F to V."""