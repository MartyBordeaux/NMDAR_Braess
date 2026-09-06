#!/usr/bin/env python3
from __future__ import annotations
import numpy as np

PULSE_COUNT = 25
PULSE_INTERVAL_MS = 5.1
GUARD_MS = 0.8
DROP_INITIAL_INTERVALS = 5

def interpulse_plateau_from_trace(time_ms, po):
    t = np.asarray(time_ms, float)
    y = np.asarray(po, float)
    if t.ndim != 1 or y.shape[-1] != len(t):
        raise ValueError("time/trace shape mismatch")
    vals = []
    for i in range(PULSE_COUNT - 1):
        if i < DROP_INITIAL_INTERVALS:
            continue
        a = i * PULSE_INTERVAL_MS + GUARD_MS
        b = (i + 1) * PULSE_INTERVAL_MS - GUARD_MS
        m = (t >= a - 1e-12) & (t <= b + 1e-12)
        if not np.any(m):
            raise RuntimeError(f"no samples in interval {i+1}: {a}-{b} ms")
        vals.append(np.nanmedian(y[..., m], axis=-1))
    interval_medians = np.stack(vals, axis=-1)
    return np.nanmedian(interval_medians, axis=-1), interval_medians

def r_plateau(time_ms, control_po, blocked_po):
    c, _ = interpulse_plateau_from_trace(time_ms, control_po)
    b, _ = interpulse_plateau_from_trace(time_ms, blocked_po)
    r = np.divide(b, c, out=np.full_like(np.asarray(b, float), np.nan),
                  where=np.asarray(c) > 1e-14)
    return c, b, r

if __name__ == "__main__":
    t = np.arange(0, 123.001, 0.025)
    c, b, r = r_plateau(t, np.full_like(t, .2), np.full_like(t, .4))
    assert abs(float(c)-.2) < 1e-12
    assert abs(float(b)-.4) < 1e-12
    assert abs(float(r)-2.0) < 1e-12
    print("plateau observable self-test PASS")
