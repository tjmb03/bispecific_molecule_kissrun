"""
The screen: find concentrations at which the tumour contact kills but does not
sustain, and normal tissue is not engaged at all.

Four layers, each answering something the one before cannot:

  binding    how much bridge forms          equilibrium, no time
  barrier    how much drug arrives          peripheral receptor sink
  kinetics   how long the contact lasts     k_off, multivalency
  thresholds what that duration produces    degranulation vs transcription

The verdict is a therapeutic window in administered concentration. Its lower
edge is where the tumour contact first holds long enough to kill; its upper
edge is where it holds long enough to drive cytokine.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import numpy as np

from .binding import trimer, penetration_factor
from .kinetics import kills, activates, _lifetime_from_trimer

__all__ = ["Candidate", "System", "screen", "sweep_koff"]


@dataclass(frozen=True)
class Candidate:
    name: str
    ka: float              # tumour arm KD, nM
    kb: float              # receptor arm KD, nM
    avidity: float         # fold tightening when both arms engage a surface
    koff_receptor: float   # receptor arm off-rate, per hour


@dataclass(frozen=True)
class System:
    antigen_tumour: float = 5.00       # nM equivalent
    antigen_normal: float = 0.25       # nM equivalent, low-level normal tissue
    receptor: float = 0.50             # nM equivalent on the engaged cell
    receptor_peripheral: float = 50.0  # nM equivalent, circulating pool (the sink)


def _longest_run(mask):
    """Start and stop indices of the longest contiguous True run.

    The window must be contiguous: a dosing range interrupted by a region of
    cytokine risk is not a therapeutic window, and reporting first-to-last
    would span the gap.
    """
    best = cur = None
    best_len = 0
    for i, v in enumerate(mask):
        if v and cur is None:
            cur = i
        elif not v and cur is not None:
            if i - cur > best_len:
                best, best_len = (cur, i - 1), i - cur
            cur = None
    if cur is not None and len(mask) - cur > best_len:
        best, best_len = (cur, len(mask) - 1), len(mask) - cur
    return best


def _trimer_profile(conc, cand: Candidate, sys: System, antigen):
    arriving = np.asarray(conc, float) * penetration_factor(
        cand.kb, sys.receptor_peripheral)
    return trimer(arriving, cand.ka, cand.kb, antigen, sys.receptor, cand.avidity)


def screen(cand: Candidate, sys: System | None = None, conc=None, **kw) -> dict:
    """Therapeutic window for one candidate.

    Three conditions must hold simultaneously:

      1. the tumour contact lasts long enough to degranulate      efficacy
      2. it does NOT last long enough to drive cytokine           safety
      3. the normal-tissue contact does not even reach killing    selectivity
    """
    sys = sys or System()
    conc = np.logspace(-5, 3, 5000) if conc is None else np.asarray(conc, float)

    t_tum = _trimer_profile(conc, cand, sys, sys.antigen_tumour)
    t_nor = _trimer_profile(conc, cand, sys, sys.antigen_normal)

    kill_tum = kills(t_tum, cand.koff_receptor, **kw)
    cyto_tum = activates(t_tum, cand.koff_receptor, **kw)
    kill_nor = kills(t_nor, cand.koff_receptor, **kw)

    ok = kill_tum & ~cyto_tum & ~kill_nor

    out = dict(asdict(cand))
    out["penetration"] = penetration_factor(cand.kb, sys.receptor_peripheral)
    out["optimum_nM"] = float(conc[np.argmax(t_tum)])
    out["kill_onset_nM"] = float(conc[np.argmax(kill_tum)]) if kill_tum.any() else float("nan")
    out["cytokine_onset_nM"] = float(conc[np.argmax(cyto_tum)]) if cyto_tum.any() else float("inf")
    out["offtumour_onset_nM"] = float(conc[np.argmax(kill_nor)]) if kill_nor.any() else float("inf")

    run = _longest_run(ok)
    if run is not None:
        lo, hi = float(conc[run[0]]), float(conc[run[1]])
        out.update(window_lo=lo, window_hi=hi,
                   margin=float(np.log10(hi / lo)), verdict="GO", reason="")
    else:
        if not kill_tum.any():
            reason = "never engages the tumour long enough to kill"
        elif (kill_tum & ~cyto_tum).sum() == 0:
            reason = "cytokine threshold reached before killing is established"
        else:
            reason = "normal tissue is engaged wherever the tumour is"
        out.update(window_lo=float("nan"), window_hi=float("nan"),
                   margin=float("nan"), verdict="NO-GO", reason=reason)
    return out


def sweep_koff(cand: Candidate, koff_values, sys: System | None = None, **kw):
    """How the therapeutic window responds to receptor-arm off-rate.

    The receptor arm sits on every side of the trade. A slower off-rate
    lengthens every contact, which brings killing within reach but also pushes
    the tumour contact past the cytokine threshold; and a tighter arm is
    sequestered peripherally and never arrives. The optimum is where those
    pressures balance.
    """
    sys = sys or System()
    rows = []
    for koff in np.asarray(koff_values, float):
        kb = cand.kb * koff / cand.koff_receptor   # k_on held fixed
        c = Candidate(cand.name, cand.ka, kb, cand.avidity, koff)
        r = screen(c, sys, **kw)
        r["halflife_min"] = 60.0 * np.log(2) / koff
        rows.append(r)
    return rows
