"""Run the screen across the panel and produce the figures.

    python scripts/run_screen.py
"""
from __future__ import annotations
import math, sys, pathlib
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from kissrun import (Candidate, System, screen, sweep_koff, DEFAULTS)
from kissrun.panel import PANEL

SYS = System()


def fmt(v, inf="unbounded", nan="—", spec=".3f"):
    if isinstance(v, float) and math.isinf(v): return inf
    if isinstance(v, float) and math.isnan(v): return nan
    return format(v, spec)


def panel_table():
    print("=" * 104)
    print("KISS-AND-RUN SELECTIVITY SCREEN")
    print(f"  kills if the contact survives {DEFAULTS['tau_kill_min']:.0f} min; "
          f"drives cytokine past {DEFAULTS['tau_cytokine_min']:.0f} min")
    print(f"  tumour antigen {SYS.antigen_tumour} nM vs normal tissue {SYS.antigen_normal} nM "
          f"({SYS.antigen_tumour/SYS.antigen_normal:.0f}-fold contrast)")
    print("=" * 104)
    hdr = (f"{'candidate':<9}{'KA':>7}{'KB':>8}{'avid':>6}{'koff/h':>8}"
           f"{'kills>=':>10}{'cyto>=':>10}{'window (nM)':>24}{'margin':>9}  verdict")
    print(hdr); print("-" * 104)
    rows = []
    for c in PANEL:
        r = screen(c, SYS); rows.append(r)
        win = (f"{fmt(r['window_lo'],spec='.4f')} - {fmt(r['window_hi'])}"
               if r["verdict"] == "GO" else "—")
        print(f"{r['name']:<9}{r['ka']:>7.1f}{r['kb']:>8.1f}{r['avidity']:>6.0f}"
              f"{r['koff_receptor']:>8.0f}{fmt(r['kill_onset_nM'],spec='.4f'):>10}"
              f"{fmt(r['cytokine_onset_nM'],inf='never',spec='.4f'):>10}{win:>24}"
              f"{fmt(r['margin'],inf='inf',spec='+.2f'):>9}  {r['verdict']} {r['reason']}")
    return rows


def koff_table():
    print("\n" + "=" * 104)
    print("RECEPTOR-ARM SWEEP — where the optimum comes from")
    print("  a slower off-rate lengthens every contact - which brings killing within")
    print("  reach, and also pushes the tumour contact past the cytokine threshold;")
    print("  and a tighter arm is sequestered peripherally and never arrives")
    print("=" * 104)
    ref = Candidate("ref", ka=16.0, kb=1.0, avidity=21.0, koff_receptor=72.0)
    koffs = [3.6, 7.2, 14.4, 36, 72, 180, 360, 720, 1800]
    print(f"{'koff /h':>10}{'contact t½':>13}{'KD (nM)':>10}{'penetration':>13}{'margin':>11}  verdict")
    print("-" * 104)
    rows = sweep_koff(ref, koffs, SYS)
    for r in rows:
        print(f"{r['koff_receptor']:>10.1f}{r['halflife_min']:>11.2f} min{r['kb']:>10.3f}"
              f"{r['penetration']:>13.4f}{fmt(r['margin'],inf='inf',spec='+.2f'):>11}  {r['verdict']}")
    return rows


def figure_mechanism():
    import subprocess, sys as _s
    subprocess.run([_s.executable, str(pathlib.Path(__file__).with_name("make_figure.py"))], check=True)
    return


if __name__ == "__main__":
    panel_table()
    koff_table()
    figure_mechanism()
