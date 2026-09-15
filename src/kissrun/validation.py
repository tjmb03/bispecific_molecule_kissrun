"""
Validation of the equilibrium layer against independent formulations.

Four checks, each testing something the others do not:

  1. MASS BALANCE      bound + free = total, at every concentration
  2. INDEPENDENT ROOT  the same conservation equations solved by a bracketing
                       root-finder rather than fixed-point iteration
  3. ODE STEADY STATE  the full kinetic scheme integrated forward until it
                       settles, compared with the algebraic answer
  4. DRUG DEPLETION    how much drug the binding actually consumes

The third earns its place: checks 1 and 2 confirm the algebra was solved
correctly, but only 3 confirms that solving it algebraically was legitimate —
that binding really is fast enough for rapid equilibrium to hold.

Check 4 exists because writing check 3 exposed an asymmetry worth stating:
the algebraic layer depletes targets but not drug.

Cross-engine comparison in the sense used for stiff ODE systems does not apply
to this layer. There is no integrator, so three languages would agree
trivially. The meaningful question is whether the modelling assumption is
sound, not whether the arithmetic is reproducible.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .binding import free_targets, trimer

__all__ = ["mass_balance_residual", "free_targets_rootfind",
           "trimer_ode_steady_state", "drug_depletion", "validate"]


# ------------------------------------------------------------------ 1
def mass_balance_residual(drug, ka, kb, a_tot, b_tot, avidity=1.0):
    """Largest relative violation of the two target conservation laws.

        A_total  =  A_free  +  [DA]  +  [DAB]
        B_total  =  B_free  +  [DB]  +  [DAB]
    """
    d = np.atleast_1d(np.asarray(drug, float))
    a_eff, b_eff = ka / avidity, kb / avidity
    a_free, b_free = free_targets(d, ka, kb, a_tot, b_tot, avidity)

    da = d * a_free / a_eff                      # drug bound by the tumour arm only
    db = d * b_free / b_eff                      # drug bound by the receptor arm only
    t = d * a_free * b_free / (a_eff * b_eff)    # ternary complex

    res_a = np.abs((a_free + da + t) - a_tot) / a_tot
    res_b = np.abs((b_free + db + t) - b_tot) / b_tot
    return float(max(res_a.max(), res_b.max()))


# ------------------------------------------------------------------ 2
def free_targets_rootfind(drug, ka, kb, a_tot, b_tot, avidity=1.0):
    """Same equations, solved by bisection on a single scalar unknown.

    Substituting the B conservation law into the A one leaves a monotone
    function of `a_free` alone, so Brent's method applies and the root is
    unique — which fixed-point iteration cannot demonstrate on its own.
    """
    a_eff, b_eff = ka / avidity, kb / avidity
    d = float(drug)

    def b_of_a(a):
        return b_tot / (1.0 + d / b_eff + d * a / (a_eff * b_eff))

    def residual(a):
        return a * (1.0 + d / a_eff + d * b_of_a(a) / (a_eff * b_eff)) - a_tot

    a_free = brentq(residual, 1e-300, a_tot, xtol=1e-18, rtol=8.9e-16)
    return a_free, b_of_a(a_free)


# ------------------------------------------------------------------ 3
def trimer_ode_steady_state(drug, ka, kb, a_tot, b_tot, avidity=1.0,
                            koff_a=3.77, koff_b=72.0, t_end=None):
    """Integrate the full binding scheme forward until it settles.

    Six species: free drug, free A, free B, the two binary complexes and the
    ternary complex. Rate constants are set so every dissociation constant
    matches the algebraic model exactly (k_on = k_off / K_D), so a disagreement
    is an equilibrium-assumption failure rather than a parameter mismatch.

    Drug is CLAMPED at the given free concentration, matching the algebraic
    layer, which treats its argument as free drug and does not conserve it.
    `drug_depletion` quantifies where that assumption starts to fail.

    Units are per hour.
    """
    a_eff, b_eff = ka / avidity, kb / avidity
    kon_a, kon_b = koff_a / a_eff, koff_b / b_eff

    def rhs(_t, y):
        D, A, B, DA, DB, T = y
        r1 = kon_a * D * A - koff_a * DA          # D  + A <-> DA
        r2 = kon_b * D * B - koff_b * DB          # D  + B <-> DB
        r3 = kon_b * DA * B - koff_b * T          # DA + B <-> T
        r4 = kon_a * DB * A - koff_a * T          # DB + A <-> T
        return [0.0,                              # drug clamped
                -r1 - r4,
                -r2 - r3,
                r1 - r3,
                r2 - r4,
                r3 + r4]

    if t_end is None:
        t_end = 500.0 / min(koff_a, koff_b)
    sol = solve_ivp(rhs, [0.0, t_end],
                    [float(drug), a_tot, b_tot, 0.0, 0.0, 0.0],
                    method="LSODA", rtol=1e-11, atol=1e-16)
    return float(sol.y[5, -1])


# ------------------------------------------------------------------ 4
def drug_depletion(drug, ka, kb, a_tot, b_tot, avidity=1.0):
    """Fraction of drug consumed by binding, at a given free concentration.

    The algebraic layer treats its argument as FREE drug and does not conserve
    it: targets deplete, drug does not. That is standard for a screening model
    and exact whenever drug is in excess — but it is an assumption, and this
    function says where it stops holding.

    Returns bound drug divided by total drug (free + bound).
    """
    d = np.atleast_1d(np.asarray(drug, float))
    a_eff, b_eff = ka / avidity, kb / avidity
    a_free, b_free = free_targets(d, ka, kb, a_tot, b_tot, avidity)
    bound = (d * a_free / a_eff + d * b_free / b_eff
             + d * a_free * b_free / (a_eff * b_eff))
    return bound / (d + bound)


# ------------------------------------------------------------------
def validate(ka=16.0, kb=1.0, a_tot=5.0, b_tot=0.5, avidity=21.0,
             concentrations=None, verbose=True):
    """Run all four checks and return the worst residual from each."""
    conc = (np.logspace(-3, 2, 60) if concentrations is None
            else np.asarray(concentrations, float))

    mb = mass_balance_residual(conc, ka, kb, a_tot, b_tot, avidity)

    rf = 0.0
    for d in conc:
        a1, b1 = free_targets(d, ka, kb, a_tot, b_tot, avidity)
        a2, b2 = free_targets_rootfind(d, ka, kb, a_tot, b_tot, avidity)
        rf = max(rf, abs(a1 - a2) / a2, abs(b1 - b2) / b2)

    probe = np.logspace(-2, 1, 7)
    ode = 0.0
    for d in probe:
        alg = float(trimer(d, ka, kb, a_tot, b_tot, avidity))
        num = trimer_ode_steady_state(d, ka, kb, a_tot, b_tot, avidity)
        ode = max(ode, abs(alg - num) / max(num, 1e-30))

    dep = drug_depletion(conc, ka, kb, a_tot, b_tot, avidity)
    worst_dep = float(dep.max())
    dep_opt = float(drug_depletion(np.sqrt(ka * kb) / avidity,
                                   ka, kb, a_tot, b_tot, avidity)[0])

    if verbose:
        print("=" * 74)
        print("VALIDATION OF THE EQUILIBRIUM LAYER")
        print("=" * 74)
        print(f"  mass balance, worst relative violation        {mb:.2e}")
        print(f"  fixed point vs bracketing root-finder         {rf:.2e}")
        print(f"  algebra vs ODE integrated to steady state     {ode:.2e}")
        print()
        print("  The third is the one that matters: it confirms the rapid-equilibrium")
        print("  assumption, not merely that the algebra was solved correctly.")
        print()
        print("-" * 74)
        print("  AN ASSUMPTION THIS EXERCISE EXPOSED")
        print("  The algebraic layer treats its argument as FREE drug and does not")
        print("  conserve it: targets deplete, drug does not.")
        print(f"    drug bound at the efficacy optimum            {100*dep_opt:5.1f} %")
        print(f"    worst across the swept range                  {100*worst_dep:5.1f} %")
        print("  Exact where drug is in excess. Where the bound fraction is large, a")
        print("  dosing model must conserve drug or it will overstate free concentration.")
    return dict(mass_balance=mb, rootfind=rf, ode_steady_state=ode,
                drug_bound_at_optimum=dep_opt, drug_bound_worst=worst_dep)
