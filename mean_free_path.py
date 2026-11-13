"""Phonon mean free path l(T) = 3 kappa / (C v), with the Debye heat capacity.

Debye parameters from the specific-heat coefficient of Y2Ti2O7 (see run_real.py);
the same lattice parameters are assumed for Tb2Ti2O7 (approximation).
"""
import numpy as np, matplotlib.pyplot as plt
from scipy.integrate import quad

KB = 1.380649e-23
v, thetaD = 3425.0, 449.0
n = 88 / (10.1e-10)**3            # atoms per m^3 (cubic cell a ~ 10.1 A, 8 f.u. x 11 atoms)

def C_debye(T):
    """Debye heat capacity per unit volume [J m^-3 K^-1]."""
    f = lambda x: x**4 / (4 * np.sinh(x / 2)**2) if x > 0 else 0.0
    return 9 * n * KB * (T / thetaD)**3 * quad(f, 0, min(thetaD / T, 50.0), limit=200)[0]

def mfp(T, k):
    return np.array([3 * kk / (C_debye(t) * v) for t, kk in zip(T, k)])

if __name__ == '__main__':
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for f, lab, c in [('data/y2ti2o7_kappa_H0.csv', 'Y$_2$Ti$_2$O$_7$', 'C0'),
                      ('data/tb2ti2o7_sampleB_kappa_H0.csv', 'Tb$_2$Ti$_2$O$_7$', 'k')]:
        d = np.genfromtxt(f, delimiter=',', names=True)
        l = mfp(d['T'], d['kappa'])
        ax.loglog(d['T'], l * 1e9, 'o', ms=3.5, c=c, label=lab)
        for T0 in (0.3, 4.0):
            print(f'{lab:22s} l({T0} K) = {np.interp(T0, d["T"], l)*1e9:9.1f} nm')
    ax.axhline(2*np.sqrt(0.63*0.15/np.pi)*1e6, ls='--', c='0.5')
    ax.text(0.35, 2*np.sqrt(0.63*0.15/np.pi)*1e6*1.3, 'Casimir length of the Y$_2$Ti$_2$O$_7$ sample', fontsize=8, color='0.4')
    ax.set(xlabel='T (K)', ylabel='phonon mean free path ℓ (nm)', title='ℓ = 3κ / (C v), Debye C')
    ax.legend(); ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout(); plt.savefig('figures/mean_free_path.png', dpi=150)
