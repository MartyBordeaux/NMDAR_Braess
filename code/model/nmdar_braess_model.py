#!/usr/bin/env python3
"""Reference implementation of the final NMDAR Braess models.

This compact publication version preserves the equations and forcing used in the
final topology-first analysis.  It is independent of the historical pipeline
layout and is intended to make the model definition directly inspectable.

Base state order (8 explicit + conserved P2):
    P0, P1A, P1B, PD1, PD2, PC1, PC2, PO
    P2 = 1 - sum(explicit states)

Mechanism index 0 = control.
Mechanism index 1 = complete B-route removal (Ro25 abstraction).

The B-slow extension uses duplicated A/B downstream states and a common
probability stock; every downstream transition of the B route is multiplied by
sB.  sB=1 nests the Base model exactly after aggregation.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
import numpy as np

PULSE_COUNT = 25
PULSE_INTERVAL_MS = 5.1
TRAIN_END_MS = (PULSE_COUNT - 1) * PULSE_INTERVAL_MS
GUARD_MS = 0.8
DROP_INITIAL_INTERVALS = 5

PARAMETER_NAMES = (
    "kon_A_per_ms", "kon_B_per_ms", "koff_A_per_ms", "koff_B_per_ms",
    "kf_plus_per_ms", "ks_plus_per_ms", "kf_minus_per_ms", "ks_minus_per_ms",
    "kd1_plus_per_ms", "kd1_minus_per_ms", "kd2_plus_per_ms", "kd2_minus_per_ms",
)


def time_grid(stop_ms: float, dt_ms: float) -> np.ndarray:
    x = np.arange(0.0, stop_ms, dt_ms, dtype=float)
    if len(x) == 0 or not math.isclose(x[-1], stop_ms, abs_tol=1e-12):
        x = np.append(x, stop_ms)
    return x


def pulse_amplitude(G_peak: float, tau_G_ms: float, G_floor: float = 0.0) -> float:
    """Per-pulse raw amplitude giving the requested final-pulse occupancy."""
    occupancy = (G_peak - G_floor) / (1.0 - G_floor)
    raw_target = occupancy / (1.0 - occupancy)
    pulse_times = np.arange(PULSE_COUNT) * PULSE_INTERVAL_MS
    accumulation = np.sum(np.exp(-(TRAIN_END_MS - pulse_times) / tau_G_ms))
    return float(raw_target / accumulation)


def forcing(t_ms: float, G_peak: float, tau_G_ms: float, G_floor: float = 0.0) -> float:
    """Normalized saturating 25-pulse forcing used by the final reruns.

    The pulse train is accumulated before 122.4 ms.  After the final pulse the
    corresponding raw accumulator decays continuously with tau_G.
    """
    amp = pulse_amplitude(G_peak, tau_G_ms, G_floor)
    occupancy = (G_peak - G_floor) / (1.0 - G_floor)
    raw_end = occupancy / (1.0 - occupancy)
    if t_ms >= TRAIN_END_MS:
        raw = raw_end * math.exp(-(t_ms - TRAIN_END_MS) / tau_G_ms)
    else:
        raw = 0.0
        for pulse in np.arange(PULSE_COUNT) * PULSE_INTERVAL_MS:
            if t_ms >= pulse - 1e-12:
                raw += amp * math.exp(-(t_ms - pulse) / tau_G_ms)
    return G_floor + (1.0 - G_floor) * raw / (1.0 + raw)


def _rate(p: dict[str, float], name: str) -> float:
    return float(p[name])


def base_rhs(y: np.ndarray, G: float, p: dict[str, float], blocked: bool = False) -> np.ndarray:
    """Derivative for the conserved nine-state Base network."""
    P0, P1A, P1B, PD1, PD2, PC1, PC2, PO = y
    P2 = 1.0 - float(np.sum(y))
    link = 0.0 if blocked else 1.0

    konA = _rate(p, "kon_A_per_ms"); konB = _rate(p, "kon_B_per_ms")
    koffA = _rate(p, "koff_A_per_ms"); koffB = _rate(p, "koff_B_per_ms")
    kfp = _rate(p, "kf_plus_per_ms"); ksp = _rate(p, "ks_plus_per_ms")
    kfm = _rate(p, "kf_minus_per_ms"); ksm = _rate(p, "ks_minus_per_ms")
    kd1p = _rate(p, "kd1_plus_per_ms"); kd1m = _rate(p, "kd1_minus_per_ms")
    kd2p = _rate(p, "kd2_plus_per_ms"); kd2m = _rate(p, "kd2_minus_per_ms")

    d = np.empty(8, dtype=float)
    d[0] = link*koffB*P1B + koffA*P1A - G*konA*P0 - link*G*konB*P0
    d[1] = G*konA*P0 - koffA*P1A + koffB*P2 - G*konB*P1A
    d[2] = link*(G*konB*P0 - koffB*P1B + koffA*P2 - G*konA*P1B)
    d[3] = kd1p*P2 - kd1m*PD1
    d[4] = kd2p*P2 - kd2m*PD2
    d[5] = kfp*P2 - kfm*PC1 + ksm*PO - ksp*PC1
    d[6] = ksp*P2 - ksm*PC2 + kfm*PO - kfp*PC2
    d[7] = ksp*PC1 - ksm*PO + kfp*PC2 - kfm*PO
    return d


def _rk4_step(rhs, y, t0, h, *args):
    g0 = args[-3](t0); gm = args[-3](t0 + 0.5*h); g1 = args[-3](t0 + h)
    fixed = args[:-3]
    k1 = rhs(y, g0, *fixed)
    k2 = rhs(y + 0.5*h*k1, gm, *fixed)
    k3 = rhs(y + 0.5*h*k2, gm, *fixed)
    k4 = rhs(y + h*k3, g1, *fixed)
    return y + h*(k1 + 2*k2 + 2*k3 + k4)/6.0


def simulate_base(p: dict[str, float], G_peak: float, tau_G_ms: float,
                  dt_ms: float = 0.1, stop_ms: float = 1500.0,
                  blocked: bool = False) -> tuple[np.ndarray, np.ndarray]:
    f = lambda t: forcing(t, G_peak, tau_G_ms)
    t = time_grid(stop_ms, dt_ms)
    y = np.zeros((len(t), 8), dtype=float); y[0, 0] = 1.0
    for i in range(len(t)-1):
        h = t[i+1] - t[i]
        g0, gm, g1 = f(t[i]), f(t[i] + 0.5*h), f(t[i+1])
        k1 = base_rhs(y[i], g0, p, blocked)
        k2 = base_rhs(y[i] + 0.5*h*k1, gm, p, blocked)
        k3 = base_rhs(y[i] + 0.5*h*k2, gm, p, blocked)
        k4 = base_rhs(y[i] + h*k3, g1, p, blocked)
        y[i+1] = y[i] + h*(k1 + 2*k2 + 2*k3 + k4)/6.0
    return t, y


def train_plateau_from_trace(t: np.ndarray, PO: np.ndarray) -> float:
    """Model analogue of the experimental negative inter-pulse plateau.

    For each adjacent pulse pair, take the median PO after 0.8-ms guards,
    discard the first five intervals, then take the median of the remaining 19.
    """
    pulse_times = np.arange(PULSE_COUNT) * PULSE_INTERVAL_MS
    vals = []
    for a, b in zip(pulse_times[:-1], pulse_times[1:]):
        m = (t >= a + GUARD_MS) & (t <= b - GUARD_MS)
        if np.any(m):
            vals.append(float(np.median(PO[m])))
    vals = vals[DROP_INITIAL_INTERVALS:]
    if len(vals) != 19:
        raise RuntimeError(f"Expected 19 retained inter-pulse intervals, got {len(vals)}")
    return float(np.median(vals))


def base_r_plateau(p: dict[str, float], G_peak: float, tau_G_ms: float,
                   dt_ms: float = 0.1) -> float:
    tc, yc = simulate_base(p, G_peak, tau_G_ms, dt_ms, stop_ms=TRAIN_END_MS, blocked=False)
    tb, yb = simulate_base(p, G_peak, tau_G_ms, dt_ms, stop_ms=TRAIN_END_MS, blocked=True)
    c = train_plateau_from_trace(tc, yc[:, 7])
    b = train_plateau_from_trace(tb, yb[:, 7])
    return b / c


# ---- B-slow extension -----------------------------------------------------
# state order: P0,P1A,P1B,P2A,P2B,PD1A,PD2A,PC1A,PC2A,POA,
#              PD1B,PD2B,PC1B,PC2B,POB


def bslow_rhs(y: np.ndarray, G: float, p: dict[str, float], sB: float,
              blocked: bool = False) -> np.ndarray:
    d = np.zeros_like(y)
    link = 0.0 if blocked else 1.0
    def edge(i, j, r):
        flux = r*y[i]; d[i] -= flux; d[j] += flux
    konA=_rate(p,"kon_A_per_ms"); konB=_rate(p,"kon_B_per_ms")
    koffA=_rate(p,"koff_A_per_ms"); koffB=_rate(p,"koff_B_per_ms")
    edge(0,1,G*konA); edge(1,0,koffA)
    edge(1,3,G*konB)
    edge(0,2,link*G*konB); edge(2,0,link*koffB)
    edge(2,4,link*G*konA)
    edge(3,1,koffB); edge(4,1,koffB)
    edge(3,2,link*koffA); edge(4,2,link*koffA)
    def downstream(p2,pd1,pd2,pc1,pc2,po,s):
        edge(p2,pd1,s*_rate(p,"kd1_plus_per_ms")); edge(pd1,p2,s*_rate(p,"kd1_minus_per_ms"))
        edge(p2,pd2,s*_rate(p,"kd2_plus_per_ms")); edge(pd2,p2,s*_rate(p,"kd2_minus_per_ms"))
        edge(p2,pc1,s*_rate(p,"kf_plus_per_ms")); edge(pc1,p2,s*_rate(p,"kf_minus_per_ms"))
        edge(p2,pc2,s*_rate(p,"ks_plus_per_ms")); edge(pc2,p2,s*_rate(p,"ks_minus_per_ms"))
        edge(pc1,po,s*_rate(p,"ks_plus_per_ms")); edge(po,pc1,s*_rate(p,"ks_minus_per_ms"))
        edge(pc2,po,s*_rate(p,"kf_plus_per_ms")); edge(po,pc2,s*_rate(p,"kf_minus_per_ms"))
    downstream(3,5,6,7,8,9,1.0)
    downstream(4,10,11,12,13,14,sB)
    return d


def simulate_bslow(p: dict[str,float], G_peak: float, tau_G_ms: float, sB: float,
                   dt_ms: float=0.1, blocked: bool=False) -> tuple[np.ndarray,np.ndarray]:
    f=lambda x: forcing(x,G_peak,tau_G_ms)
    t=time_grid(TRAIN_END_MS,dt_ms)
    y=np.zeros((len(t),15),float); y[0,0]=1.0
    for i in range(len(t)-1):
        h=t[i+1]-t[i]; g0=f(t[i]); gm=f(t[i]+0.5*h); g1=f(t[i+1])
        k1=bslow_rhs(y[i],g0,p,sB,blocked)
        k2=bslow_rhs(y[i]+0.5*h*k1,gm,p,sB,blocked)
        k3=bslow_rhs(y[i]+0.5*h*k2,gm,p,sB,blocked)
        k4=bslow_rhs(y[i]+h*k3,g1,p,sB,blocked)
        y[i+1]=y[i]+h*(k1+2*k2+2*k3+k4)/6.0
    return t,y


def bslow_r_plateau(p: dict[str,float], G_peak: float, tau_G_ms: float, sB: float,
                    dt_ms: float=0.1) -> float:
    tc,yc=simulate_bslow(p,G_peak,tau_G_ms,sB,dt_ms,False)
    tb,yb=simulate_bslow(p,G_peak,tau_G_ms,sB,dt_ms,True)
    c=train_plateau_from_trace(tc,yc[:,9]+yc[:,14])
    b=train_plateau_from_trace(tb,yb[:,9]+yb[:,14])
    return b/c


if __name__ == "__main__":
    print("NMDAR Braess reference model; import this module to run simulations.")
