"""Regenerate figures/mechanism.png with a schematic of the mechanism itself."""
from __future__ import annotations
import sys, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from kissrun import (Candidate, System, screen, sweep_koff, trimer,
                     penetration_factor, contact_lifetime, DEFAULTS)
from kissrun.panel import PANEL

TEAL, ROSE, AMBER = "#0E6E73", "#A8324A", "#B4690E"
GREY, INK, MUTE, RULE = "#9AA7AE", "#10212B", "#62727C", "#DCE3E6"
WASH, BLUSH = "#EEF4F4", "#FBF1F3"
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.edgecolor": RULE,
                     "axes.linewidth": 0.9, "xtick.color": MUTE,
                     "ytick.color": MUTE, "xtick.labelsize": 9,
                     "ytick.labelsize": 9})
OUT = pathlib.Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(exist_ok=True)
SYS = System()

fig = plt.figure(figsize=(13.2, 7.4))
gs = fig.add_gridspec(2, 3, height_ratios=[0.92, 1.0], hspace=0.40, wspace=0.28,
                      left=0.055, right=0.985, top=0.885, bottom=0.075)

# ============================ A — the schematic ============================
ax = fig.add_subplot(gs[0, :]); ax.axis("off")
ax.set_xlim(0, 20); ax.set_ylim(-0.1, 6.2)
ax.set_title("A · Two dwell times, one contact", fontsize=12, color=INK,
             fontweight="bold", loc="left", pad=4)

def cell(x, y, r, fc, ec, label, sub=None):
    ax.add_patch(Circle((x, y), r, facecolor=fc, edgecolor=ec, lw=1.8, zorder=2))
    ax.text(x, y + 0.08, label, fontsize=9, color=ec, ha="center", va="center",
            fontweight="bold", zorder=3)
    if sub:
        ax.text(x, y + r + 0.34, sub, fontsize=8.2, color=MUTE, ha="center")

def bridges(x0, x1, y, n, col):
    if n == 1:
        ys = [y]
    else:
        ys = np.linspace(y - 0.62, y + 0.62, n)
    for yy in ys:
        ax.plot([x0, x1], [yy, yy], color=col, lw=2.0, solid_capstyle="round", zorder=1)
        ax.plot([(x0 + x1) / 2], [yy], marker="o", ms=4.5, color=col, zorder=3)

# --- left: low density -> kiss and run
ax.add_patch(plt.Rectangle((0.3, 0.55), 9.0, 4.7, facecolor=WASH,
                           edgecolor=TEAL, lw=1.0, zorder=0, alpha=0.75))
cell(2.4, 3.1, 1.15, WASH, TEAL, "T cell")
cell(7.0, 3.1, 1.15, "#FFFFFF", TEAL, "tumour", "brief engagement")
bridges(3.55, 5.85, 3.1, 2, TEAL)
ax.text(4.7, 4.02, "2 bridges", fontsize=9, color=MUTE, ha="center")

ax.text(4.8, 1.20, "3.3 min  —  degranulates, then moves on",
        fontsize=10, color=TEAL, ha="center", fontweight="bold")
ax.text(4.8, 0.78, "past τ$_{kill}$, short of τ$_{cytokine}$  —  the therapeutic window", fontsize=8.6,
        color=MUTE, ha="center")

# --- right: high density -> holds and signals
ax.add_patch(plt.Rectangle((10.5, 0.55), 9.2, 4.7, facecolor=BLUSH,
                           edgecolor=ROSE, lw=1.0, zorder=0, alpha=0.55))
cell(12.6, 3.1, 1.15, WASH, TEAL, "T cell")
cell(17.3, 3.1, 1.15, "#FFFFFF", ROSE, "tumour", "sustained engagement")
bridges(13.75, 16.15, 3.1, 5, ROSE)
ax.text(14.95, 4.02, "5 bridges", fontsize=9, color=MUTE, ha="center")
ax.text(14.95, 1.20, "3.6 h  —  transcription, and cytokine",
        fontsize=10, color=ROSE, ha="center", fontweight="bold")
ax.text(14.95, 0.78, "past τ$_{cytokine}$  —  sustained activation",
        fontsize=8.6, color=MUTE, ha="center")

ax.text(10.0, 5.75, "granules are pre-formed; cytokine needs transcription  —  so two thresholds, not one",
        fontsize=9.5, color=INK, ha="center", style="italic")
ax.text(10.0, 0.18, "illustrative counts, p$_{unbound}$ = 0.25, receptor k$_{off}$ = 72/h",
        fontsize=8.2, color=MUTE, ha="center", style="italic")

# ============================ B — exponential ==============================
ax = fig.add_subplot(gs[1, 0])
n = np.arange(1, 9)
for p, c in ((0.15, TEAL), (0.25, AMBER), (0.35, ROSE)):
    ax.semilogy(n, contact_lifetime(n, 72.0, p_unbound=p), "o-", color=c,
                lw=2.2, ms=5, label=f"p$_{{unbound}}$ = {p}")
ax.axhline(DEFAULTS["tau_kill_min"], color=TEAL, lw=1.2, ls=(0, (4, 3)))
ax.axhline(DEFAULTS["tau_cytokine_min"], color=ROSE, lw=1.2, ls=(0, (4, 3)))
ax.text(7.9, DEFAULTS["tau_kill_min"] * 1.5, "τ$_{kill}$", fontsize=9, color=TEAL, ha="right")
ax.text(7.9, DEFAULTS["tau_cytokine_min"] * 1.5, "τ$_{cytokine}$", fontsize=9, color=ROSE, ha="right")
ax.set_xlabel("simultaneous bridges", fontsize=10, color=INK)
ax.set_ylabel("contact lifetime (min, log)", fontsize=10, color=INK)
ax.set_title("B · Lifetime is exponential in bridge count", fontsize=11.5,
             color=INK, fontweight="bold", loc="left", pad=9)
ax.legend(frameon=False, fontsize=8.5, loc="upper left")
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(True, color="#EEF2F4", lw=0.8); ax.set_axisbelow(True)

# ============================ C — lead candidate ===========================
ax = fig.add_subplot(gs[1, 1])
lead = PANEL[0]
conc = np.logspace(-4, 3, 2000)
pen = penetration_factor(lead.kb, SYS.receptor_peripheral)
for dens, lab, col in ((SYS.antigen_tumour, "tumour", TEAL),
                       (SYS.antigen_normal, "normal tissue", ROSE)):
    t = trimer(conc * pen, lead.ka, lead.kb, dens, SYS.receptor, lead.avidity)
    ax.loglog(conc, np.clip(contact_lifetime(t * DEFAULTS["synapse_gain"],
                                             lead.koff_receptor), 1e-3, 1e12),
              color=col, lw=2.6, label=lab)
ax.axhline(DEFAULTS["tau_kill_min"], color=TEAL, lw=1.2, ls=(0, (4, 3)))
ax.axhline(DEFAULTS["tau_cytokine_min"], color=ROSE, lw=1.2, ls=(0, (4, 3)))
r = screen(lead, SYS)
if r["verdict"] == "GO":
    hi = r["window_hi"] if np.isfinite(r["window_hi"]) else conc.max()
    ax.axvspan(r["window_lo"], hi, color=TEAL, alpha=0.10, lw=0)
    ax.text(np.sqrt(r["window_lo"] * hi), 2e8, "therapeutic\nwindow",
            fontsize=8.8, color=TEAL, ha="center", linespacing=1.4)
ax.set_xlabel("administered concentration (nM, log)", fontsize=10, color=INK)
ax.set_ylabel("contact lifetime (min, log)", fontsize=10, color=INK)
ax.set_title(f"C · {lead.name} — kill without sustaining", fontsize=11.5,
             color=INK, fontweight="bold", loc="left", pad=9)
ax.legend(frameon=False, fontsize=9, loc="lower right")
ax.set_ylim(1e-2, 1e10)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(True, color="#EEF2F4", lw=0.8); ax.set_axisbelow(True)

# ============================ D — the barrier ==============================
ax = fig.add_subplot(gs[1, 2])
ref = Candidate("ref", ka=16.0, kb=1.0, avidity=21.0, koff_receptor=72.0)
koffs = np.logspace(np.log10(2.0), np.log10(3000), 45)
rows = sweep_koff(ref, koffs, SYS)
ax.semilogx([r["koff_receptor"] for r in rows],
            [r["penetration"] for r in rows], color=AMBER, lw=2.6)
fails = [r["koff_receptor"] for r in rows if r["verdict"] == "NO-GO"]
if fails:
    ax.axvspan(koffs.min(), max(fails), color=ROSE, alpha=0.11, lw=0)
    ax.text(max(fails) * 0.82, 0.62, "too tight —\nnever arrives", fontsize=9,
            color=ROSE, ha="right", linespacing=1.4)
ax.text(2200, 0.30, "too loose —\ncontact too brief", fontsize=9, color=MUTE,
        ha="right", linespacing=1.4)
ax.set_xlabel("receptor-arm k$_{off}$ (1/h, log)", fontsize=10, color=INK)
ax.set_ylabel("fraction reaching the tumour", fontsize=10, color=INK)
ax.set_title("D · Why affinity has an optimum", fontsize=11.5, color=INK,
             fontweight="bold", loc="left", pad=9)
ax.set_ylim(0, 1.02)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(True, color="#EEF2F4", lw=0.8); ax.set_axisbelow(True)

fig.suptitle("Kinetic selectivity: dwell time, not occupancy", x=0.055,
             ha="left", fontsize=13, color=INK, fontweight="bold", y=0.965)
fig.text(0.055, 0.928,
         "Killing needs a brief contact; cytokine needs a sustained one. The window is the gap between them.",
         fontsize=10, color=MUTE)

fig.savefig(OUT / "mechanism.png", dpi=200, facecolor="white")
print("wrote", OUT / "mechanism.png")
