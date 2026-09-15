"""
Equilibrium ternary-complex formation for a T-cell-engaging bispecific.

Rapid-binding equilibrium with target depletion, solved by fixed-point
iteration. No time variable: this layer answers "how much bridge forms at
this concentration", not "how long does it last".

The kinetic layer in kinetics.py supplies the second question.
"""
from __future__ import annotations
import numpy as np

__all__ = ["free_targets", "trimer", "peak_concentration", "penetration_factor"]

_MAX_ITER = 500
_TOL = 1e-14


def free_targets(drug, ka, kb, a_tot, b_tot, avidity=1.0):
    """Free tumour antigen and free receptor at equilibrium.

    Parameters
    ----------
    drug : float or array
        Free bispecific concentration (nM).
    ka, kb : float
        Dissociation constants of the tumour arm and the receptor arm (nM).
    a_tot, b_tot : float
        Total accessible tumour antigen and receptor (nM equivalent).
    avidity : float
        Fold tightening of both effective affinities once the assembled
        complex engages a surface with both arms. 1.0 means no avidity.

    Returns
    -------
    (a_free, b_free) : same shape as `drug`
    """
    d = np.asarray(drug, dtype=float)
    a_eff, b_eff = ka / avidity, kb / avidity
    a_free = np.full_like(d, a_tot, dtype=float)
    b_free = np.full_like(d, b_tot, dtype=float)
    for _ in range(_MAX_ITER):
        a_new = a_tot / (1.0 + d / a_eff + d * b_free / (a_eff * b_eff))
        b_new = b_tot / (1.0 + d / b_eff + d * a_free / (a_eff * b_eff))
        if np.max(np.abs(a_new - a_free)) < _TOL and np.max(np.abs(b_new - b_free)) < _TOL:
            a_free, b_free = a_new, b_new
            break
        a_free, b_free = a_new, b_new
    return a_free, b_free


def trimer(drug, ka, kb, a_tot, b_tot, avidity=1.0):
    """Productive ternary complex concentration (nM)."""
    d = np.asarray(drug, dtype=float)
    a_eff, b_eff = ka / avidity, kb / avidity
    a_free, b_free = free_targets(d, ka, kb, a_tot, b_tot, avidity)
    return d * a_free * b_free / (a_eff * b_eff)


def peak_concentration(ka, kb, avidity=1.0):
    """Closed-form optimum: the geometric mean of the effective arm affinities.

    Exact in the target-excess limit. Depends only on the PRODUCT of the two
    constants, so an asymmetric design peaks where a balanced one does.
    Use a numerical maximum of `trimer` when depletion is material.
    """
    return np.sqrt(ka * kb) / avidity


def penetration_factor(kb, receptor_peripheral):
    """Fraction of drug surviving peripheral receptor binding en route to tumour.

    The binding-site barrier: a bispecific with a very tight receptor arm binds
    the first T cell it meets and never reaches the tumour core. Modelled as a
    single-site sink on the receptor arm, which has no avidity benefit because
    nothing is bridged in transit.

    This term is what prevents the screen from recommending unbounded affinity.
    """
    return 1.0 / (1.0 + receptor_peripheral / kb)
