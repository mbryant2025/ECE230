from constants import *
from sympy.physics.units import cm, volt, second
import sympy as sp

def calc_mu_nN(Nd):
    """Calculate the electron mobility in a n-type region

    Parameters
    ----------
    Nd : float
        The donor concentration in cm^-3.

    Returns
    -------
    mu_nN : float
        The electron mobility in the n-type region in cm^2/Vs.
    """

    mu_nN = 92 + 1268 / (1 + (remove_units(Nd) / 1.3e17)**(0.91))

    return mu_nN * cm**2 / volt / second

def calc_mu_pN(Nd):
    """Calculate the hole mobility in a p-type region

    Parameters
    ----------
    Nd : float
        The acceptor concentration in cm^-3.

    Returns
    -------
    mu_pN : float
        The hole mobility in the p-type region in cm^2/Vs.
    """

    mu_pN = 130 + 370 / (1 + (remove_units(Nd) / 8.0e17)**(1.25))

    return mu_pN * cm**2 / volt / second

def calc_mu_nP(Na):
    """Calculate the electron mobility in a p-type region

    Parameters
    ----------
    Na : float
        The acceptor concentration in cm^-3.

    Returns
    -------
    mu_nP : float
        The electron mobility in the n-type region in cm^2/Vs.
    """

    mu_nP = 232 + 1180 / (1 + (remove_units(Na) / 8.0e16)**(0.90))

    return mu_nP * cm**2 / volt / second

def calc_mu_pP(Na):
    """Calculate the hole mobility in a p-type region

    Parameters
    ----------
    Na : float
        The acceptor concentration in cm^-3.

    Returns
    -------
    mu_pP : float
        The hole mobility in the p-type region in cm^2/Vs.
    """

    mu_pP = 49.7 + 418.3 / (1 + (remove_units(Na) / 1.6e17)**(0.70))

    return mu_pP * cm**2 / volt / second

def calc_tau_rec_P(Na):
    """Calculate the recombination time in a p-type region

    Parameters
    ----------
    Na : float
        The acceptor concentration in cm^-3.

    Returns
    -------
    tau_rec_P : float
        The recombination time in the p-type region in seconds.
    """

    Na = remove_units(Na)

    tau_rec_P = 1 / (3.45e-12 * Na+ 9.50e-32 * Na**2)

    return tau_rec_P * second

def calc_tau_rec_N(Nd):
    """Calculate the recombination time in a n-type region

    Parameters
    ----------
    Nd : float
        The donor concentration in cm^-3.

    Returns
    -------
    tau_rec_N : float
        The recombination time in the n-type region in seconds.
    """

    Nd = remove_units(Nd)

    tau_rec_N = 1 / (7.80e-13 * Nd+ 1.80e-31 * Nd**2)

    return tau_rec_N * second

def calc_tau_gen_N(Nd):
    """Calculate the generation time in a n-type region

    Parameters
    ----------
    Nd : float
        The donor concentration in cm^-3.

    Returns
    -------
    tau_gen_N : float
        The generation time in the n-type region in seconds.
    """

    tau_rec_N = calc_tau_rec_N(Nd)

    return tau_rec_N / 75

def calc_tau_gen_P(Na):
    """Calculate the generation time in a p-type region

    Parameters
    ----------
    Na : float
        The acceptor concentration in cm^-3.

    Returns
    -------
    tau_gen_P : float
        The generation time in the p-type region in seconds.
    """

    tau_rec_P = calc_tau_rec_P(Na)

    return tau_rec_P / 75

def calc_L_nP(Na):
    """Calculate the diffusion length of minority electrons

    Parameters
    ----------
    Na: float
        The acceptor concentration in cm^-3.

    Returns
    -------
    L_np : float
        The diffusion length of minority electrons in the p-type region in cm.
    """

    mu_nP = calc_mu_nP(Na)
    D_nP = calc_D_nP(mu_nP)
    tau_nP = calc_tau_gen_P(Na)

    return sp.sqrt(D_nP * tau_nP)

def calc_L_pN(Nd):
    """Calculate the diffusion length of minority holes

    Parameters
    ----------
    Nd: float
        The donor concentration in cm^-3.

    Returns
    -------
    L_pn : float
        The diffusion length of minority holes in the n-type region in cm.
    """

    mu_pN = calc_mu_pN(Nd)
    D_pN = calc_D_pN(mu_pN)
    tau_pN = calc_tau_gen_N(Nd)

    return sp.sqrt(D_pN * tau_pN)

def calc_D_pN(mu_pN):
    """Calculate the hole diffusion coefficient in the n-type region

    Parameters
    ----------
    mu_pN : float
        The hole mobility in the n-type region in cm^2/Vs.

    Returns
    -------
    D_pN : float
        The hole diffusion coefficient in the n-type region in cm^2/s.
    """

    return mu_pN * k_BT_300K_DIVIDED_q

def calc_D_nP(mu_nP):
    """Calculate the electron diffusion coefficient in the p-type region

    Parameters
    ----------
    mu_nP : float
        The electron mobility in the p-type region in cm^2/Vs.

    Returns
    -------
    D_nP : float
        The electron diffusion coefficient in the p-type region in cm^2/s.
    """

    return mu_nP * k_BT_300K_DIVIDED_q

