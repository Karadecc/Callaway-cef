"""
Fit kappa(T[, H]) data with the Callaway + CEF model.

Workflow used in the README:
  stage 1  fit the non-magnetic analogue (e.g. Y2Ti2O7)  -> L, A, B  (no resonant term)
  stage 2  fix the lattice terms, fit the Tb compound at H = 0 -> C_a_b (+ level energies, bounded)
  stage 3  predict kappa(T, H) with NO new free parameter and compare with field data
           (then, optionally, a joint fit over all fields)

CLI:
  python fit.py configs/tb2ti2o7_template.json --stage reference --out lattice.json
  python fit.py configs/tb2ti2o7_template.json --stage magnetic --lattice lattice.json --out cef.json
"""
from __future__ import annotations

import argparse
import json

import numpy as np
from scipy.optimize import least_squares

from callaway_cef import kappa_data, scheme_from_config


def fit_kappa(T, k, H, p0: dict, free: list, scheme=None, bounds: dict | None = None):
    """Least squares on log(kappa) (data spans decades) with parameters in log10.

    Returns best params, 1-sigma relative errors, correlation matrix, result object.
    """
    bounds = bounds or {}
    names = list(free)
    x0 = np.log10([p0[n] for n in names])
    lo = np.array([np.log10(bounds[n][0]) if n in bounds else -np.inf for n in names])
    hi = np.array([np.log10(bounds[n][1]) if n in bounds else np.inf for n in names])
    x0 = np.clip(x0, lo + 1e-9, hi - 1e-9)

    def unpack(x):
        p = dict(p0)
        p.update({n: 10.0**xi for n, xi in zip(names, x)})
        return p

    def resid(x):
        return np.log(kappa_data(T, H, unpack(x), scheme)) - np.log(k)

    res = least_squares(resid, x0, bounds=(lo, hi), x_scale="jac")
    dof = max(len(k) - len(names), 1)
    s2 = 2 * res.cost / dof
    cov = np.linalg.pinv(res.jac.T @ res.jac) * s2           # in log10 space
    sig = np.sqrt(np.clip(np.diag(cov), 0, None))
    corr = cov / np.outer(sig, sig).clip(1e-300)
    rel_err = {n: np.log(10) * s for n, s in zip(names, sig)}  # d p / p
    return unpack(res.x), rel_err, corr, res


def report(names, p, rel_err, corr, res, n):
    rms = np.sqrt(2 * res.cost / n)
    print(f"  rms log-residual = {rms:.3f}  (~{100*rms:.1f}% typical deviation)")
    for nm in names:
        print(f"  {nm:>8s} = {p[nm]:.4g}  +- {100*rel_err[nm]:.1f}%")
    if len(names) > 1:
        print("  correlation matrix (|r| > 0.9 means the data cannot separate them):")
        print("           " + " ".join(f"{nm:>8s}" for nm in names))
        for nm, row in zip(names, corr):
            print(f"  {nm:>8s} " + " ".join(f"{r:8.2f}" for r in row))


def load_csv(path):
    """CSV with header; columns T [K], kappa [W/m/K], optional H [T]."""
    d = np.genfromtxt(path, delimiter=",", names=True)
    H = d["H"] if "H" in d.dtype.names else np.zeros_like(d["T"])
    return d["T"], d["kappa"], H


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--stage", choices=["reference", "magnetic"], required=True)
    ap.add_argument("--out", default=None, help="write fitted params (json)")
    ap.add_argument("--lattice", default=None, help="json from the reference stage (L, A, B)")
    a = ap.parse_args()
    cfg = json.load(open(a.config))
    st = cfg[a.stage]
    lattice = json.load(open(a.lattice)) if a.lattice else {}
    p0 = {**cfg["fixed"], **lattice, **st["init"]}
    scheme = scheme_from_config(cfg) if a.stage == "magnetic" else None
    T, k, H = load_csv(st["data"])
    if "Hmax" in st:
        keep = H <= st["Hmax"]
        T, k, H = T[keep], k[keep], H[keep]
    p, err, corr, res = fit_kappa(T, k, H, p0, st["free"], scheme, st.get("bounds"))
    print(f"[{a.stage}] {st['data']}  ({len(T)} points)")
    report(st["free"], p, err, corr, res, len(T))
    if a.out:
        json.dump({n: p[n] for n in st["free"]}, open(a.out, "w"), indent=2)


if __name__ == "__main__":
    main()
