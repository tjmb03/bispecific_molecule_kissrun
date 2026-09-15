"""Tests for kissrun.

The analytical assertions are the point: they hold regardless of parameter
values, so a change that breaks them is a change to the physics, not a
tuning difference.
"""
import math
import numpy as np
import pytest

from kissrun import (trimer, free_targets, peak_concentration, penetration_factor,
                     bonds_per_contact, contact_lifetime, threshold_bonds,
                     kills, activates, DEFAULTS,
                     Candidate, System, screen, sweep_koff)
from kissrun.panel import PANEL


# ---------------------------------------------------------------- binding

def test_trimer_peaks_at_geometric_mean_in_target_excess():
    """The closed form: peak at sqrt(KA*KB), exact when target is plentiful."""
    ka, kb = 16.0, 1.0
    d = np.logspace(-3, 3, 20000)
    # target excess: densities far above the drug concentrations of interest
    t = trimer(d, ka, kb, a_tot=1e5, b_tot=1e5)
    assert math.isclose(d[np.argmax(t)], peak_concentration(ka, kb),
                        rel_tol=2e-3)


def test_peak_depends_only_on_the_product_of_arm_affinities():
    """Arms of 0.3 and 30 peak where arms of 3 and 3 do."""
    d = np.logspace(-3, 3, 20000)
    asym = trimer(d, 0.3, 30.0, 1e5, 1e5)
    symm = trimer(d, 3.0, 3.0, 1e5, 1e5)
    assert math.isclose(d[np.argmax(asym)], d[np.argmax(symm)], rel_tol=5e-3)


def test_density_moves_peak_height_not_position():
    """Second analytical result: two separable levers."""
    d = np.logspace(-3, 3, 20000)
    low = trimer(d, 16.0, 1.0, 1e5, 1e3)
    high = trimer(d, 16.0, 1.0, 1e5, 1e4)
    assert high.max() > low.max()
    assert math.isclose(d[np.argmax(low)], d[np.argmax(high)], rel_tol=5e-3)


def test_avidity_shifts_the_optimum_by_its_own_factor():
    assert math.isclose(peak_concentration(16.0, 1.0, avidity=21.0),
                        peak_concentration(16.0, 1.0) / 21.0, rel_tol=1e-12)


def test_free_target_never_exceeds_total():
    a, b = free_targets(np.logspace(-3, 3, 200), 16.0, 1.0, 5.0, 0.5, 21.0)
    assert np.all(a <= 5.0 + 1e-12) and np.all(a > 0)
    assert np.all(b <= 0.5 + 1e-12) and np.all(b > 0)


def test_depletion_is_material_at_the_optimum():
    """The closed form assumes target excess. Here it does not hold."""
    d = peak_concentration(16.0, 1.0, 21.0)
    _, b_free = free_targets(d, 16.0, 1.0, 5.0, 0.5, 21.0)
    assert b_free / 0.5 < 0.5, "expected substantial receptor sequestration"


# ---------------------------------------------------------------- barrier

def test_penetration_falls_as_the_receptor_arm_tightens():
    tight = penetration_factor(kb=0.1, receptor_peripheral=50.0)
    loose = penetration_factor(kb=100.0, receptor_peripheral=50.0)
    assert tight < loose < 1.0


# ---------------------------------------------------------------- kinetics

def test_single_bond_lifetime_is_the_reciprocal_off_rate():
    assert math.isclose(contact_lifetime(1, koff_per_hr=72.0), 60.0 / 72.0)


def test_lifetime_is_exponential_in_bond_count():
    """Each extra bond multiplies lifetime by 1/p, not adds to it."""
    p = 0.25
    lives = [contact_lifetime(n, 72.0, p_unbound=p) for n in (1, 2, 3, 4)]
    ratios = [b / a for a, b in zip(lives, lives[1:])]
    assert all(math.isclose(r, 1 / p, rel_tol=1e-9) for r in ratios)


def test_threshold_bonds_inverts_the_lifetime_relation():
    koff, tau = 72.0, 5.0
    n = threshold_bonds(koff, tau)
    assert math.isclose(contact_lifetime(n, koff), tau, rel_tol=1e-9)


def test_killing_threshold_is_reached_before_the_cytokine_one():
    """The whole basis of the window: degranulation needs less dwell than
    transcription, so there is a range where a contact kills and moves on."""
    assert DEFAULTS["tau_kill_min"] < DEFAULTS["tau_cytokine_min"]
    n_kill = threshold_bonds(72.0, DEFAULTS["tau_kill_min"])
    n_cyto = threshold_bonds(72.0, DEFAULTS["tau_cytokine_min"])
    assert n_kill < n_cyto


def test_a_contact_between_the_thresholds_kills_without_activating():
    """Pick a trimer level landing between the two dwell times."""
    koff = 72.0
    n_mid = 0.5 * (threshold_bonds(koff, DEFAULTS["tau_kill_min"])
                   + threshold_bonds(koff, DEFAULTS["tau_cytokine_min"]))
    t_mid = n_mid / DEFAULTS["synapse_gain"]
    assert bool(kills(t_mid, koff))
    assert not bool(activates(t_mid, koff))


def test_bonds_scale_linearly_with_trimer():
    assert math.isclose(bonds_per_contact(0.2) / bonds_per_contact(0.1), 2.0)


# ---------------------------------------------------------------- screen

def test_every_panel_candidate_returns_a_verdict():
    for c in PANEL:
        r = screen(c)
        assert r["verdict"] in ("GO", "NO-GO")
        if r["verdict"] == "GO":
            assert r["window_lo"] > 0
            assert r["window_hi"] > r["window_lo"]


def test_losing_avidity_pushes_the_window_to_higher_concentrations():
    """BSP-01 and BSP-02 are the same binder; only the geometry differs.

    Avidity lengthens every contact, so it brings killing within reach at much
    lower concentration - and, in the two-threshold picture, also brings the
    cytokine threshold closer. The window moves, it does not simply widen.
    """
    lead = screen(PANEL[0])
    no_avidity = screen(PANEL[1])
    assert no_avidity["window_lo"] > lead["window_lo"]


def test_the_window_is_bounded_by_the_two_thresholds():
    """Lower edge is where killing starts, upper edge where cytokine does."""
    r = screen(PANEL[0])
    assert r["verdict"] == "GO"
    assert math.isclose(r["window_lo"], r["kill_onset_nM"], rel_tol=1e-6)
    assert r["window_hi"] < r["cytokine_onset_nM"]


def test_the_window_is_contiguous():
    """A range interrupted by cytokine risk is not a therapeutic window."""
    from kissrun.kinetics import kills as _k, activates as _a
    from kissrun.binding import trimer as _t, penetration_factor as _p
    c, s = PANEL[0], System()
    r = screen(c, s)
    if r["verdict"] == "GO":
        probe = np.logspace(np.log10(r["window_lo"]), np.log10(r["window_hi"]), 200)
        arriving = probe * _p(c.kb, s.receptor_peripheral)
        t = _t(arriving, c.ka, c.kb, s.antigen_tumour, s.receptor, c.avidity)
        assert np.all(_k(t, c.koff_receptor))
        assert not np.any(_a(t, c.koff_receptor))


def test_an_over_tight_receptor_arm_fails():
    """The binding-site barrier is what stops the screen preferring infinite affinity."""
    c = Candidate("too-tight", ka=16.0, kb=0.02, avidity=21.0, koff_receptor=1.44)
    r = screen(c)
    assert r["verdict"] == "NO-GO"
    assert r["reason"], "a NO-GO must say why"
    assert r["penetration"] < 0.01, "an over-tight arm should barely reach the tumour"


def test_every_go_reports_both_threshold_crossings():
    for c in PANEL:
        r = screen(c)
        if r["verdict"] == "GO":
            assert r["kill_onset_nM"] < r["cytokine_onset_nM"]


def test_koff_sweep_has_both_failure_modes():
    """Too tight fails on delivery; the sweep must not be monotone in affinity."""
    c = Candidate("ref", ka=16.0, kb=1.0, avidity=21.0, koff_receptor=72.0)
    rows = sweep_koff(c, [3.6, 7.2, 14.4, 72.0, 360.0, 1800.0])
    verdicts = [r["verdict"] for r in rows]
    assert "NO-GO" in verdicts and "GO" in verdicts


def test_system_defaults_are_a_twenty_fold_density_contrast():
    s = System()
    assert math.isclose(s.antigen_tumour / s.antigen_normal, 20.0)


# ---------------------------------------------------------------- validation

from kissrun import (mass_balance_residual, free_targets_rootfind,
                     trimer_ode_steady_state, drug_depletion, validate)

_SYS = dict(ka=16.0, kb=1.0, a_tot=5.0, b_tot=0.5, avidity=21.0)


def test_mass_balance_holds_to_machine_precision():
    conc = np.logspace(-3, 2, 200)
    assert mass_balance_residual(conc, **_SYS) < 1e-12


def test_fixed_point_agrees_with_an_independent_root_finder():
    """Confirms the iteration converged to the right root, not merely to a root."""
    for d in np.logspace(-3, 2, 25):
        a1, b1 = free_targets(d, _SYS["ka"], _SYS["kb"], _SYS["a_tot"],
                              _SYS["b_tot"], _SYS["avidity"])
        a2, b2 = free_targets_rootfind(d, **_SYS)
        assert math.isclose(float(a1), a2, rel_tol=1e-10)
        assert math.isclose(float(b1), b2, rel_tol=1e-10)


def test_algebra_matches_the_kinetic_scheme_at_steady_state():
    """The check that validates the modelling assumption, not just the arithmetic.

    Integrating the full six-species binding scheme forward must land on the
    algebraic answer. If it does not, rapid equilibrium is not a safe
    simplification for these rate constants.
    """
    for d in np.logspace(-2, 1, 5):
        alg = float(trimer(d, _SYS["ka"], _SYS["kb"], _SYS["a_tot"],
                           _SYS["b_tot"], _SYS["avidity"]))
        num = trimer_ode_steady_state(d, **_SYS)
        assert math.isclose(alg, num, rel_tol=1e-8)


def test_drug_depletion_is_reported_and_material():
    """The algebra does not conserve drug. This asserts the size of that gap.

    At the efficacy optimum most of the drug is bound, so the free-drug
    assumption is badly violated there. A dosing model must conserve drug.
    """
    opt = peak_concentration(_SYS["ka"], _SYS["kb"], _SYS["avidity"])
    bound = float(drug_depletion(opt, **_SYS)[0])
    assert 0.5 < bound < 1.0, "expected substantial drug consumption at the optimum"


def test_drug_depletion_vanishes_when_drug_is_in_excess():
    """The assumption is exact in the limit it was made for."""
    tiny_target = dict(_SYS, a_tot=1e-6, b_tot=1e-7)
    bound = float(drug_depletion(1.0, **tiny_target)[0])
    assert bound < 1e-4


def test_validate_returns_all_residuals_below_tolerance():
    r = validate(verbose=False)
    assert r["mass_balance"] < 1e-12
    assert r["rootfind"] < 1e-10
    assert r["ode_steady_state"] < 1e-8


# ---------------------------------------------------------------- identifiability

from kissrun import trimer_thresholds, equivalent_parameterisation, report


def test_kinetic_layer_reduces_to_two_trimer_thresholds():
    """Five kinetic parameters, two identifiable combinations."""
    t_kill, t_cyt = trimer_thresholds(72.0)
    assert 0 < t_kill < t_cyt


def test_an_equivalent_parameterisation_gives_identical_thresholds():
    """Constructive demonstration of the degeneracy, not an assertion about it."""
    alt = equivalent_parameterisation(72.0, synapse_gain_new=80.0, p_unbound_new=0.5)
    a = trimer_thresholds(72.0)
    b = trimer_thresholds(72.0, alt["tau_kill_min"], alt["tau_cytokine_min"],
                          alt["p_unbound"], alt["synapse_gain"])
    assert math.isclose(a[0], b[0], rel_tol=1e-12)
    assert math.isclose(a[1], b[1], rel_tol=1e-12)


def test_degenerate_parameters_produce_bit_identical_windows():
    """The degeneracy is not approximate: every verdict and edge matches."""
    alt = equivalent_parameterisation(72.0, synapse_gain_new=80.0, p_unbound_new=0.5)
    for c in PANEL:
        a = screen(c)
        b = screen(c, **alt)
        assert a["verdict"] == b["verdict"]
        if a["verdict"] == "GO":
            assert math.isclose(a["window_lo"], b["window_lo"], rel_tol=1e-12)
            assert math.isclose(a["window_hi"], b["window_hi"], rel_tol=1e-12)
            assert math.isclose(a["margin"], b["margin"], rel_tol=1e-12)


def test_report_returns_a_vanishing_difference():
    r = report(verbose=False)
    assert r["max_relative_difference"] < 1e-12
