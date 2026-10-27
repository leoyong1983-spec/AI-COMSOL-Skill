"""Shared finite-table and scale-aware numerical comparisons."""
import math
from reaction_core import InputError


def verify_readback(before, after, relative_tolerance=1e-8, absolute_tolerance=1e-14):
    """Compare each SI-valued column on its own magnitude, including trace fluxes."""
    if not before or len(before)!=len(after):raise InputError("save/reopen row count mismatch")
    worst=0.0
    legacy=0.0
    for row, reread in zip(before, after):
        if not row or len(row)!=len(reread):raise InputError("save/reopen table shape mismatch")
        for a,b in zip(row,reread):
            if not math.isfinite(a) or not math.isfinite(b):raise InputError("nonfinite save/reopen value")
            delta=abs(a-b)
            worst=max(worst,delta/(absolute_tolerance+relative_tolerance*max(abs(a),abs(b))))
            legacy=max(legacy,delta/max(1,abs(a)))
    if worst>1:raise InputError("save/reopen mismatch on a column's own magnitude; tolerance ratio="+str(worst))
    return {"save_reopen_max_scaled_difference":legacy,"save_reopen_max_tolerance_ratio":worst,
            "save_reopen_relative_tolerance":relative_tolerance,"save_reopen_absolute_tolerance_SI":absolute_tolerance}
