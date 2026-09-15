"""A synthetic candidate panel, built to exercise the levers.

Not real molecules. Each differs from BSP-01 in one property so the screen's
response to that property is visible.
"""
from .screen import Candidate

PANEL = [
    Candidate("BSP-01", ka=16.0, kb=1.0,  avidity=21.0, koff_receptor=72.0),
    Candidate("BSP-02", ka=16.0, kb=1.0,  avidity=1.0,  koff_receptor=72.0),
    Candidate("BSP-03", ka=5.0,  kb=0.2,  avidity=8.0,  koff_receptor=14.4),
    Candidate("BSP-04", ka=0.2,  kb=5.0,  avidity=8.0,  koff_receptor=360.0),
    Candidate("BSP-05", ka=1.0,  kb=1.0,  avidity=6.0,  koff_receptor=72.0),
    Candidate("BSP-06", ka=16.0, kb=100.0,avidity=21.0, koff_receptor=720.0),
]
