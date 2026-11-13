"""
Callaway model of phonon thermal conductivity with resonant scattering on
crystal-electric-field (CEF) levels, and its magnetic-field dependence.

    kappa(T) = (kB / 2 pi^2 v) (kB T / hbar)^3  Int_0^{thetaD/T} tau_c(x) x^4 e^x / (e^x - 1)^2 dx

    1/tau_c = v/L                                   boundaries (Casimir)
            + D w                                   dislocations (strain field), optional
            + A w^4                                 point defects (Rayleigh)
            + B w^2 T exp(-thetaD / (b T))          umklapp
            + sum_ij C_ij w^4 / [(w_ij^2 - w^2)^2 + (g w_ij w)^2] * F_ij(T, H)
                                                    resonant scattering on CEF level pairs

All energies (phonon energy w = hbar omega / kB, level energies) are in kelvin,
so D is in s^-1 K^-1, A in s^-1 K^-4, B in s^-1 K^-3, C_ij in s^-1. Fields in tesla.
N-processes are not included (Callaway's N-term is usually negligible for
strongly scattering, glass-like systems; see README for when this matters).
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field

import numpy as np

KB = 1.380649e-23          # J/K
HBAR = 1.054571817e-34     # J s
MUB_OVER_KB = 0.67171381563  # mu_B / kB in K/T

_trapz = getattr(np, "trapezoid", None) or np.trapz


@dataclass
class Multiplet:
    """A CEF singlet (deg=1) or (quasi-)doublet (deg=2) of the rare-earth ion.

    E0    : centre energy above the ground multiplet at H = 0 [K]
    delta : zero-field splitting of a non-Kramers quasi-doublet [K]
    m     : effective moment along the local easy axis [mu_B];
            a doublet splits as  E0 +- sqrt((delta/2)^2 + (m muB H cos(theta))^2)
    """
    E0: float
    deg: int = 2
    delta: float = 0.0
    m: float = 0.0

    def energies(self, H: float, proj: float = 1.0) -> np.ndarray:
        if self.deg == 1:
            return np.array([self.E0])
        half = np.hypot(0.5 * self.delta, self.m * proj * MUB_OVER_KB * H)
        return np.array([self.E0 - half, self.E0 + half])


@dataclass
class LevelScheme:
    """CEF level scheme + which multiplet pairs scatter phonons resonantly.

    transitions : list of (a, b) multiplet indices. (a, a) is the transition
                  inside a doublet; (a, b) couples every sublevel of a to every
                  sublevel of b with the same strength C_a_b (simplification:
                  real matrix elements differ between sublevels).
    sites       : list of (weight, cos(theta)) for magnetically inequivalent
                  sites, theta = angle between H and the local easy axis.
                  Pyrochlore, H || [111]: [(0.25, 1.0), (0.75, 1/3)]
                  Pyrochlore, H || [001]: [(1.0, 1/sqrt(3))]
                  Pyrochlore, H || [1-10]: [(0.5, sqrt(2/3)), (0.5, 0.0)]
    population  : "difference" -> F = p_lower - p_upper  (net resonant absorption)
                  "sum"        -> F = p_lower + p_upper  (both levels scatter)
                  The literature uses both; which one fits is a physics question.
    width       : fractional linewidth g of each resonance (regularises w_ij -> 0).
    """
    multiplets: list
    transitions: list
    sites: list = field(default_factory=lambda: [(1.0, 1.0)])
    population: str = "difference"
    width: float = 0.05

    def with_params(self, p: dict) -> "LevelScheme":
        """Return a copy where keys E{k}, delta{k}, m{k} in p override level k."""
        s = copy.deepcopy(self)
        for k, mlt in enumerate(s.multiplets):
            for attr in ("E0", "delta", "m"):
                key = f"{'E' if attr == 'E0' else attr}{k}"
                if key in p:
                    setattr(mlt, attr, p[key])
        return s

    def _sublevels(self, H, proj):
        E, owner = [], []
        for k, mlt in enumerate(self.multiplets):
            e = mlt.energies(H, proj)
            E.extend(e)
            owner.extend([k] * len(e))
        return np.array(E), np.array(owner)

    def resonances(self, T, H, proj):
        """List of (w_ij [K], F_ij(T) array, transition key)."""
        E, owner = self._sublevels(H, proj)
        E = E - E.min()
        boltz = np.exp(-E[None, :] / T[:, None])
        p = boltz / boltz.sum(axis=1, keepdims=True)
        out = []
        for a, b in self.transitions:
            key = f"C_{a}_{b}"
            ia = np.where(owner == a)[0]
            ib = np.where(owner == b)[0]
            pairs = ([(ia[0], ia[1])] if a == b and len(ia) == 2
                     else [(i, j) for i in ia for j in ib if i != j])
            for i, j in pairs:
                lo, hi = (i, j) if E[i] <= E[j] else (j, i)
                F = p[:, lo] - p[:, hi] if self.population == "difference" else p[:, lo] + p[:, hi]
                out.append((E[hi] - E[lo], F, key))
        return out

    def resonant_rate(self, w, T, H, p):
        """Site-averaged resonant scattering rate [s^-1], shape of w (nT, nx)."""
        rate = np.zeros_like(w)
        for weight, proj in self.sites:
            for w0, F, key in self.resonances(T, H, proj):
                C = p.get(key, 0.0)
                if C == 0.0:
                    continue
                den = (w0**2 - w**2) ** 2 + (self.width * w0 * w) ** 2 + 1e-12
                rate += weight * C * w**4 / den * F[:, None]
        return rate


def kappa(T, H, p, scheme: LevelScheme | None = None, nx: int = 4000, xcut: float = 50.0,
          return_rates: bool = False):
    """Phonon thermal conductivity [W m^-1 K^-1] at temperatures T [K], field H [T].

    p : dict with v [m/s], thetaD [K], L [m], A, B, b, and C_a_b for resonances
        (plus optional E{k}, delta{k}, m{k} overriding the level scheme).
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    xmax = np.minimum(p["thetaD"] / T, xcut)
    u = np.linspace(0.0, 1.0, nx)[1:]
    x = xmax[:, None] * u[None, :]
    w = x * T[:, None]                                  # phonon energy in K

    rates = {
        "boundary": np.full_like(w, p["v"] / p["L"]),
        "dislocation": p.get("D", 0.0) * w,
        "point": p.get("A", 0.0) * w**4,
        "umklapp": p.get("B", 0.0) * w**2 * T[:, None] * np.exp(-p["thetaD"] / (p.get("b", 2.0) * T[:, None])),
    }
    if scheme is not None:
        s = scheme.with_params(p)
        rates["resonant"] = s.resonant_rate(w, T, H, p)
    total = sum(rates.values())

    weight = x**4 / (4.0 * np.sinh(x / 2.0) ** 2)       # x^4 e^x / (e^x-1)^2, overflow-safe
    integral = _trapz(weight / total, x, axis=1)
    k = KB / (2 * np.pi**2 * p["v"]) * (KB * T / HBAR) ** 3 * integral
    return (k, x, w, rates) if return_rates else k


def kappa_data(T, H, p, scheme=None):
    """kappa for data where each point may have its own field H (arrays, same length)."""
    T, H = np.asarray(T, float), np.broadcast_to(np.asarray(H, float), np.shape(T))
    out = np.empty_like(T)
    for h in np.unique(H):
        sel = H == h
        out[sel] = kappa(T[sel], h, p, scheme)
    return out


def scheme_from_config(cfg: dict) -> LevelScheme:
    lv = cfg["levels"]
    return LevelScheme(
        multiplets=[Multiplet(**{k: v for k, v in m.items() if not k.startswith("_")}) for m in lv["multiplets"]],
        transitions=[tuple(t) for t in lv["transitions"]],
        sites=[tuple(s) for s in lv.get("sites", [(1.0, 1.0)])],
        population=lv.get("population", "difference"),
        width=lv.get("width", 0.05),
    )
