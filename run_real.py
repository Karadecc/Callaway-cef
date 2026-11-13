"""Three-stage analysis on the digitised Li et al. (2013) data."""
import numpy as np, json
from callaway_cef import kappa, kappa_data, LevelScheme, Multiplet
from fit import fit_kappa, report

casimir = lambda w, t: 2*np.sqrt(w*t/np.pi)*1e-3       # mm -> m
L_Y = casimir(0.63, 0.15)                                  # Y2Ti2O7 sample
L_A = casimir(0.66, 0.16)                                  # Tb2Ti2O7 sample A
DEBYE = dict(v=3425.0, thetaD=449.0)                       # from beta of Y2Ti2O7 (Li et al., ref. 28)

# ---- stage 1: Y2Ti2O7
d = np.genfromtxt('data/y2ti2o7_kappa_H0.csv', delimiter=',', names=True)
p0 = dict(**DEBYE, L=L_Y, D=1e7, A=1e2, B=1e5, b=4.0)
free1 = ['D', 'A', 'B', 'b']
lat, e1, c1, r1 = fit_kappa(d['T'], d['kappa'], np.zeros_like(d['T']), p0, free1)
print('STAGE 1  Y2Ti2O7, L fixed to Casimir length %.2f mm' % (L_Y*1e3)); report(free1, lat, e1, c1, r1, len(d))

# ---- stage 2: Tb2Ti2O7 sample A, H = 0
def scheme(pop='difference', m0=6.0, m1=3.0, sites=((0.25, 1.0), (0.75, 1/3))):
    return LevelScheme([Multiplet(0.0, 2, 1.0, m0), Multiplet(18.0, 2, 0.0, m1)],
                       [(0, 0), (0, 1)], sites=list(sites), population=pop)

f = np.genfromtxt('data/tb2ti2o7_kappa.csv', delimiter=',', names=True)
z = f['H'] == 0
pTb0 = {**lat, 'L': L_A, 'C_0_0': 1e9, 'C_0_1': 1e10, 'E1': 18.0, 'delta0': 1.0}
free2 = ['C_0_0', 'C_0_1', 'E1', 'delta0']
bounds2 = {'E1': (15.0, 21.0), 'delta0': (0.05, 4.0)}
fits = {}
for pop in ('difference', 'sum'):
    pTb, e2, c2, r2 = fit_kappa(f['T'][z], f['kappa'][z], f['H'][z], pTb0, free2, scheme(pop), bounds2)
    print(f'\nSTAGE 2  Tb2Ti2O7 (A), H = 0, population = {pop}'); report(free2, pTb, e2, c2, r2, z.sum())
    fits[pop] = pTb

json.dump({'lattice': lat, 'tb': fits}, open('results_real.json', 'w'), indent=1, default=float)
