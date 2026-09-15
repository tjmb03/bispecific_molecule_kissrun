"""kissrun - a kinetic therapeutic-window screen for T-cell-engaging bispecifics.

Four layers: equilibrium binding, a peripheral binding-site barrier, a
multivalent contact-lifetime calculation, and two dwell-time thresholds -
degranulation and cytokine transcription.

The output is a therapeutic window in concentration. Its lower edge is where
the tumour contact first holds long enough to kill; its upper edge is where it
holds long enough to drive cytokine. That is kiss-and-run: kill and move on.
"""
from .binding import trimer, free_targets, peak_concentration, penetration_factor
from .kinetics import (bonds_per_contact, contact_lifetime, kills, activates,
                       threshold_bonds, DEFAULTS)
from .screen import Candidate, System, screen, sweep_koff
from .identifiability import trimer_thresholds, equivalent_parameterisation, report
from .validation import (mass_balance_residual, free_targets_rootfind,
                         trimer_ode_steady_state, drug_depletion, validate)

__version__ = "0.1.0"
__all__ = ["trimer","free_targets","peak_concentration","penetration_factor",
           "bonds_per_contact","contact_lifetime","kills","activates","threshold_bonds",
           "DEFAULTS","Candidate","System","screen","sweep_koff",
           "mass_balance_residual","free_targets_rootfind",
           "trimer_ode_steady_state","drug_depletion","validate",
           "trimer_thresholds","equivalent_parameterisation","report"]
