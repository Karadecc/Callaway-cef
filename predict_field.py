"""
Stage 3: predict kappa(T, H) from the zero-field fit, with no new free parameter,
and overlay the measured field data.

  python predict_field.py configs/tb2ti2o7_template.json \
         --lattice lattice.json --cef cef.json --out figures/field_prediction.png
"""
import argparse
import json

import matplotlib.pyplot as plt
import numpy as np

from callaway_cef import kappa, kappa_data, scheme_from_config
from fit import load_csv

ap = argparse.ArgumentParser()
ap.add_argument("config")
ap.add_argument("--lattice", required=True)
ap.add_argument("--cef", required=True)
ap.add_argument("--out", default="figures/field_prediction.png")
a = ap.parse_args()

cfg = json.load(open(a.config))
p = {**cfg["fixed"], **json.load(open(a.lattice)), **json.load(open(a.cef))}
scheme = scheme_from_config(cfg)
T, k, H = load_csv(cfg["magnetic"]["data"])

fields = np.unique(H)
r = np.log(kappa_data(T, H, p, scheme)) - np.log(k)
for h in fields:
    sel = H == h
    print(f"H = {h:5.2f} T : rms log-residual = {np.sqrt(np.mean(r[sel]**2)):.3f}"
          f"{'   (fitted)' if h == 0 else '   (predicted)'}")

Tf = np.logspace(np.log10(T.min()), np.log10(T.max()), 200)
fig, ax = plt.subplots(figsize=(6, 4.5))
for h, c in zip(fields, plt.cm.viridis(np.linspace(0, 0.9, len(fields)))):
    sel = H == h
    ax.loglog(T[sel], k[sel], "o", ms=4, c=c)
    ax.loglog(Tf, kappa(Tf, h, p, scheme), c=c, label=f"{h:g} T")
ax.set(xlabel="T (K)", ylabel="κ (W m$^{-1}$ K$^{-1}$)",
       title="Points: data. Lines: H = 0 fit extrapolated in field")
ax.legend(fontsize=8)
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig(a.out, dpi=150, bbox_inches="tight")
