"""
What this model can and cannot have calibrated.

The kinetic layer looks like five parameters - synapse gain, unbound
probability, receptor off-rate, and the two dwell thresholds. It is not. Every
screen verdict depends on them only through TWO derived quantities.

    the condition          g * T  >=  n*(tau, p, k_off)
    rearranged             T      >=  n* / g

so the whole kinetic machinery collapses into a threshold on trimer:

    T_kill      =  n*(tau_kill,     p, k_off) / g
    T_cytokine  =  n*(tau_cytokine, p, k_off) / g

FIVE parameters, TWO identifiable combinations, THREE degrees of freedom that
no screen outcome can recover. Raise the synapse gain and lengthen the dwell
thresholds to compensate, and every window is bit-for-bit identical.

WHY THIS IS WORTH KNOWING RATHER THAN HIDING
--------------------------------------------
It says which experiments would change an answer and which would not.

Measuring p_unbound precisely, or pinning down the synapse gain, does NOT
improve the screen: both enter only through the two thresholds. What WOULD
improve it is calibrating T_kill and T_cytokine directly against a molecule
with known clinical behaviour - two numbers, not five.

It also bounds the claim. The kinetic construction is a mechanistic account of
WHY two thresholds exist and why they differ. It is not evidence that the
particular values of g, p or tau used here are right, and nothing in the screen
output could make it so.
"""
from __future__ import annotations
import numpy as np

from .kinetics import DEFAULTS, threshold_bonds

__all__ = ["trimer_thresholds", "equivalent_parameterisation", "report"]


def trimer_thresholds(koff_per_hr, tau_kill_min=None, tau_cytokine_min=None,
                      p_unbound=None, synapse_gain=None):
    """The two quantities the screen actually depends on.

    Returns (T_kill, T_cytokine) in nM - the trimer concentrations at which a
    contact first survives each dwell threshold.
    """
    g = DEFAULTS["synapse_gain"] if synapse_gain is None else synapse_gain
    tk = DEFAULTS["tau_kill_min"] if tau_kill_min is None else tau_kill_min
    tc = (DEFAULTS["tau_cytokine_min"] if tau_cytokine_min is None
          else tau_cytokine_min)
    p = DEFAULTS["p_unbound"] if p_unbound is None else p_unbound
    return (float(threshold_bonds(koff_per_hr, tk, p) / g),
            float(threshold_bonds(koff_per_hr, tc, p) / g))


def equivalent_parameterisation(koff_per_hr, synapse_gain_new, p_unbound_new,
                                **base):
    """Construct a different (g, p, tau) triple giving identical screen output.

    Pick any new gain and unbound probability; this returns the two dwell
    thresholds that keep both trimer thresholds unchanged. Demonstrates the
    degeneracy constructively rather than asserting it.
    """
    t_kill, t_cyt = trimer_thresholds(koff_per_hr, **base)
    out = {}
    for name, t_thresh in (("tau_kill_min", t_kill),
                           ("tau_cytokine_min", t_cyt)):
        n_needed = t_thresh * synapse_gain_new
        out[name] = float(60.0 / (koff_per_hr * p_unbound_new ** (n_needed - 1.0)))
    out["synapse_gain"] = float(synapse_gain_new)
    out["p_unbound"] = float(p_unbound_new)
    return out


def report(koff_per_hr=72.0, verbose=True):
    """State the identifiable combinations and demonstrate the degeneracy."""
    t_kill, t_cyt = trimer_thresholds(koff_per_hr)
    alt = equivalent_parameterisation(koff_per_hr, 80.0, 0.5)
    t_kill2, t_cyt2 = trimer_thresholds(
        koff_per_hr, alt["tau_kill_min"], alt["tau_cytokine_min"],
        alt["p_unbound"], alt["synapse_gain"])

    if verbose:
        print("=" * 74)
        print("IDENTIFIABILITY OF THE KINETIC LAYER")
        print("=" * 74)
        print("  Five parameters enter: synapse gain, unbound probability,")
        print("  receptor off-rate, and two dwell thresholds.")
        print("  Every screen verdict depends on them only through two numbers.\n")
        print(f"{'':>26}{'T_kill (nM)':>14}{'T_cytokine (nM)':>18}")
        print(f"{'as configured':>26}{t_kill:>14.7f}{t_cyt:>18.7f}")
        print(f"{'an alternative':>26}{t_kill2:>14.7f}{t_cyt2:>18.7f}")
        print(f"\n  alternative: g = {alt['synapse_gain']:.0f}, "
              f"p = {alt['p_unbound']:.2f}, "
              f"tau_kill = {alt['tau_kill_min']:.3f} min, "
              f"tau_cytokine = {alt['tau_cytokine_min']:.3f} min")
        print("  -> different parameters, identical thresholds, identical windows.\n")
        print("  WHAT FOLLOWS")
        print("  Measuring p_unbound or the synapse gain would not change any verdict.")
        print("  Calibrating T_kill and T_cytokine against a molecule with known")
        print("  clinical behaviour would change every one. Two numbers, not five.")
    return dict(T_kill=t_kill, T_cytokine=t_cyt, alternative=alt,
                max_relative_difference=max(abs(t_kill - t_kill2) / t_kill,
                                            abs(t_cyt - t_cyt2) / t_cyt))
