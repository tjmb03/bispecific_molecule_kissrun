"""
Contact lifetime, and the two dwell times that separate killing from cytokine.

A cell-cell contact held by several independent bridges ruptures only when all
of them are momentarily unbound at the same instant:

    k_rupture = k_off * p_unbound**(n - 1)

so lifetime is EXPONENTIAL in bridge count while bridge count is only LINEAR
in target density.

WHY TWO THRESHOLDS
------------------
Cytolytic granule release and cytokine production do not require the same
engagement time. Granules are pre-formed and pre-loaded, so degranulation
follows quickly once signalling starts. Mass cytokine production requires a
transcriptional programme, which needs a sustained signal.

That gives two dwell times on the same axis:

    tau_kill      short   enough contact to trigger degranulation
    tau_cytokine  long    enough contact to drive the transcriptional programme

A contact between them kills and moves on - serial killing without sustained
activation. That is kiss-and-run, and it is the rationale for the deliberately
weak receptor arms used in current designs.

An equilibrium model cannot represent any of this. It has no k_off, so every
complex counts the same regardless of how long it survives.
"""
from __future__ import annotations
import numpy as np

__all__ = ["bonds_per_contact", "contact_lifetime", "kills", "activates",
           "threshold_bonds", "DEFAULTS"]

DEFAULTS = dict(
    synapse_gain=40.0,       # bridges per nM of trimer concentrated at the interface
    p_unbound=0.25,          # momentary unbound probability of a single bridge
    tau_kill_min=2.0,        # contact needed to trigger degranulation
    tau_cytokine_min=20.0,   # contact needed to drive cytokine transcription
)


def bonds_per_contact(trimer_nM, synapse_gain=None):
    """Simultaneous bridges at one cell-cell interface.

    The gain converts a bulk trimer concentration into a count at the
    interface. It absorbs synapse area, receptor mobility and the local
    concentrating effect, and is the least well determined input here.
    """
    g = DEFAULTS["synapse_gain"] if synapse_gain is None else synapse_gain
    return np.asarray(trimer_nM, dtype=float) * g


def contact_lifetime(n_bonds, koff_per_hr, p_unbound=None):
    """Mean contact duration in minutes.

    Ruptures only when every bridge is simultaneously off. At n = 1 this
    reduces to the single-bond lifetime 1 / k_off.
    """
    p = DEFAULTS["p_unbound"] if p_unbound is None else p_unbound
    n = np.maximum(np.asarray(n_bonds, dtype=float), 1.0)
    rate = koff_per_hr * np.power(p, n - 1.0)
    with np.errstate(divide="ignore"):
        return 60.0 / rate


def _lifetime_from_trimer(trimer_nM, koff_per_hr, **kw):
    n = bonds_per_contact(trimer_nM, kw.get("synapse_gain"))
    return contact_lifetime(n, koff_per_hr, kw.get("p_unbound"))


def kills(trimer_nM, koff_per_hr, tau_kill_min=None, **kw):
    """Does the contact last long enough to degranulate?"""
    tau = DEFAULTS["tau_kill_min"] if tau_kill_min is None else tau_kill_min
    return _lifetime_from_trimer(trimer_nM, koff_per_hr, **kw) >= tau


def activates(trimer_nM, koff_per_hr, tau_cytokine_min=None, **kw):
    """Does the contact last long enough to drive cytokine transcription?"""
    tau = (DEFAULTS["tau_cytokine_min"] if tau_cytokine_min is None
           else tau_cytokine_min)
    return _lifetime_from_trimer(trimer_nM, koff_per_hr, **kw) >= tau


def threshold_bonds(koff_per_hr, tau_min=None, p_unbound=None):
    """Bridges required for a contact to survive tau_min.

    Solves  60 / (k_off * p**(n-1))  =  tau  for n.
    """
    tau = DEFAULTS["tau_kill_min"] if tau_min is None else tau_min
    p = DEFAULTS["p_unbound"] if p_unbound is None else p_unbound
    return 1.0 + np.log(60.0 / (tau * koff_per_hr)) / np.log(p)
