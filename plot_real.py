"""Figures for the analysis of the digitised Li et al. (2013) data (run run_real.py first)."""
import numpy as np, json, matplotlib.pyplot as plt
from callaway_cef import kappa, LevelScheme, Multiplet
R = json.load(open('results_real.json')); lat = R['lattice']; fits = R['tb']
L_A = 2*np.sqrt(0.66*0.16/np.pi)*1e-3
def scheme(pop, m0=6.0, m1=3.0):
    return LevelScheme([Multiplet(0, 2, 1.0, m0), Multiplet(18, 2, 0, m1)], [(0, 0), (0, 1)],
                       sites=[(0.25, 1), (0.75, 1/3)], population=pop)
y = np.genfromtxt('data/y2ti2o7_kappa_H0.csv', delimiter=',', names=True)
f = np.genfromtxt('data/tb2ti2o7_kappa.csv', delimiter=',', names=True)
rr = np.genfromtxt('data/tb2ti2o7_ratio_vs_H.csv', delimiter=',', names=True)
st = dict(difference=dict(c='C3', ls='-', lab='F = p$_l$ − p$_u$'), sum=dict(c='C1', ls='--', lab='F = p$_l$ + p$_u$'))

fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
Tf = np.logspace(np.log10(0.3), np.log10(87), 200)
a = ax[0]
a.loglog(y['T'], y['kappa'], 'o', mfc='none', c='C0', ms=4, label='Y$_2$Ti$_2$O$_7$ (data)')
a.loglog(Tf, kappa(Tf, 0, lat), c='C0', label='Callaway fit (L fixed)')
z = f['H'] == 0
a.loglog(f['T'][z], f['kappa'][z], 's', c='k', ms=4, label='Tb$_2$Ti$_2$O$_7$ A, 0 T (data)')
Tt = np.logspace(np.log10(0.3), np.log10(25), 200)
for pop, p in fits.items():
    a.loglog(Tt, kappa(Tt, 0, {**lat, **p}, scheme(pop)), color=st[pop]['c'], ls=st[pop]['ls'], label='CEF fit, ' + st[pop]['lab'])
a.set(xlabel='T (K)', ylabel='κ (W m$^{-1}$ K$^{-1}$)', title='Stages 1–2: zero-field fits'); a.legend(fontsize=7)

a = ax[1]
Hs = np.linspace(0, 14, 141)
s = rr['T'] == 0.36
a.plot(rr['H'][s], rr['ratio'][s], 'ks', ms=5, label='data, 0.36 K')
for pop, p in fits.items():
    P = {**lat, **p}; k0 = kappa(0.36, 0, P, scheme(pop))[0]
    a.plot(Hs, [kappa(0.36, h, P, scheme(pop))[0]/k0 for h in Hs], color=st[pop]['c'], ls=st[pop]['ls'], label='prediction, ' + st[pop]['lab'])
a.set(xlabel='μ$_0$H || [111] (T)', ylabel='κ(H)/κ(0)', title='Stage 3: 0.36 K, no new parameter'); a.legend(fontsize=7)

a = ax[2]
P = {**lat, **fits['sum']}
for T0, c in zip([0.36, 0.7, 1.4, 5.0], plt.cm.plasma(np.linspace(0, 0.8, 4))):
    s = rr['T'] == T0
    a.plot(rr['H'][s], rr['ratio'][s], 'o', c=c, ms=4)
    k0 = kappa(T0, 0, P, scheme('sum'))[0]
    a.plot(Hs, [kappa(T0, h, P, scheme('sum'))[0]/k0 for h in Hs], c=c, label=f'{T0} K')
a.set_yscale('log'); a.set(xlabel='μ$_0$H || [111] (T)', ylabel='κ(H)/κ(0)',
       title='Points: data. Lines: prediction (sum)'); a.legend(fontsize=7)
for a in ax: a.spines[['top', 'right']].set_visible(False)
plt.tight_layout(); plt.savefig('figures/real_data_analysis.png', dpi=150, bbox_inches='tight')
