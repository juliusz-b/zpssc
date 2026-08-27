"""s26_laws_concept.py - concept panels for the closed-form laws.

  (a) Law A mechanism: a grating read through an upstream neighbour's
      transmission notch; the notch is asymmetric across the line, so the
      fitted peak is pulled by delta. Reflectivity drawn exaggerated.
  (b) The placement map that Law A implies: |delta(Delta)| with the
      tolerance line, the forbidden interval [Delta_lo, Delta_hi], and two
      operating windows of a +/-10 K pair: co-tuned (straddles the
      forbidden zone) and far-detuned (clear of it).
  (c) Law B mechanism: the mean multiple-access background acts as a pure
      stretch of the wavelength axis; arrows grow linearly with position,
      and the two axis references pin the stretch away.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import figstyle as FS

FS.apply()

SIG = 106.2          # pm
C_A = (4.0 / 3.0) * np.sqrt(2.0 / 3.0)

fig, ax = plt.subplots(1, 3, figsize=(7.1, 2.25))

# ------------------------- (a) Law A mechanism ------------------------------
a = ax[0]
nu = np.linspace(-420, 420, 900)
Rx = 0.35                                  # exaggerated for visibility
DET = 130.0                                # neighbour detuned by Delta*
g = np.exp(-0.5 * (nu / SIG) ** 2)
T = (1.0 - Rx * np.exp(-0.5 * ((nu - DET) / SIG) ** 2)) ** 2
seen = g * T
a.fill_between(nu / 1000, g, color='0.85', lw=0)
a.plot(nu / 1000, g, color='0.55', lw=0.8, label='true line')
a.plot(nu / 1000, T, color=FS.BLUE, lw=0.9, ls=(0, (4, 2)),
       label='upstream notch')
a.plot(nu / 1000, seen, color=FS.VERM, lw=1.2, label='line as read')
mu_seen = nu[np.argmax(seen)] / 1000
a.axvline(0, color='0.4', lw=0.6, ls=(0, (2, 2)))
a.axvline(mu_seen, color=FS.VERM, lw=0.6, ls=(0, (2, 2)))
a.annotate('', xy=(mu_seen, 1.04), xytext=(0, 1.04),
           arrowprops=dict(arrowstyle='-|>', lw=0.9, color='0.15'))
a.text(mu_seen / 2, 1.07, r'$\delta$', fontsize=8, ha='center', color='0.15')
a.text(DET / 1000, 0.30, r'neighbour at $\Delta$', fontsize=6.4,
       color=FS.BLUE, ha='center')
a.set_xlim(-0.42, 0.42); a.set_ylim(0, 1.16)
a.set_xlabel(r'$\lambda$ offset [nm]')
a.set_ylabel('reflectance, transmission [norm.]')
a.set_title('(a) Law A: the notch pulls the fit', fontsize=9)
a.legend(fontsize=6.2, loc='upper left')

# ------------------------- (b) the placement map ----------------------------
b = ax[1]
R = 0.10
eps = 1.0
D = np.linspace(0, 700, 800)
bias = C_A * R * D * np.exp(-D ** 2 / (3 * SIG ** 2))
b.plot(D, bias, color=FS.VERM, lw=1.2)
b.axhline(eps, color='0.3', lw=0.7, ls=(0, (4, 2)))
b.text(660, eps * 1.15, r'tolerance $\epsilon$', fontsize=6.4, ha='right',
       color='0.3')
# roots
from scipy.optimize import brentq
f = lambda x: C_A * R * x * np.exp(-x ** 2 / (3 * SIG ** 2)) - eps
dlo = brentq(f, 1, SIG * np.sqrt(1.5))
dhi = brentq(f, SIG * np.sqrt(1.5), 700)
b.axvspan(dlo, dhi, color=FS.VERM, alpha=0.10, lw=0)
b.text((dlo + dhi) / 2, 7.4, 'forbidden', fontsize=6.6, ha='center',
       color=FS.VERM)
for dv, lab in [(dlo, r'$\Delta_\mathrm{lo}$'), (dhi, r'$\Delta_\mathrm{hi}$')]:
    b.axvline(dv, color=FS.VERM, lw=0.6, ls=(0, (2, 2)))
    b.text(dv, -0.85, lab, fontsize=7, ha='center', va='top',
           color=FS.VERM)
# operating windows of a +/-10 K pair
S = 100.0
b.plot([0, 2 * S], [5.6, 5.6], color='0.25', lw=3, solid_capstyle='butt')
b.text(S + 45, 6.15, r'co-tuned pair, $D=0$', fontsize=6.2, ha='center',
       color='0.25')
b.plot([550 - 2 * S, 550 + 2 * S], [4.2, 4.2], color=FS.GREEN, lw=3,
       solid_capstyle='butt')
b.text(550, 4.55, r'far-detuned, $D=550$ pm', fontsize=6.2, ha='center',
       color=FS.GREEN)
b.set_xlim(0, 700); b.set_ylim(0, 9.4)
b.set_xlabel(r'pair detuning $\Delta$ [pm]')
b.set_ylabel(r'pairwise bias $|\delta|$ [pm]')
b.set_title('(b) the windows a range sweeps', fontsize=9)

# ------------------------- (c) Law B mechanism ------------------------------
c = ax[2]
W = 200.0
pos = np.linspace(-W, W, 9)
K, N = 32, 127
stretch = (K - 1) / N
for p in pos:
    c.plot([p, p], [0, 0.55], color='0.75', lw=1.2)
    dx = stretch * p              # true scale: (K-1)/N per unit position
    if abs(dx) > 2:
        c.annotate('', xy=(p + dx, 0.72), xytext=(p, 0.72),
                   arrowprops=dict(arrowstyle='-|>', lw=0.8,
                                   color=FS.VERM))
refs = [-150, 150]
for r0 in refs:
    c.plot([r0, r0], [0, 0.55], color=FS.BLUE, lw=2.0)
    c.plot([r0], [0.02], marker='^', ms=5, color=FS.BLUE)
c.text(0, 0.90, r'mean MAI bias $=$ stretch by $1+(K{-}1)/N$',
       fontsize=6.6, ha='center', color=FS.VERM)
c.text(0, -0.30, 'two pinned references remove the stretch\n'
       r'(9.4 $\rightarrow$ 2.3 pm at $K=32$, $N=127$)',
       fontsize=6.4, ha='center', color=FS.BLUE)
c.set_xlim(-230, 230); c.set_ylim(-0.55, 1.05)
c.set_xticks([-200, 0, 200])
c.set_yticks([])
c.set_xlabel(r'position in the band $\nu_k$ [pm]')
c.set_title('(c) Law B: an odd pull, pinned by refs', fontsize=9)
for sp in ('left',):
    c.spines[sp].set_visible(False)

plt.tight_layout()
fig.savefig('figs/fig_s26_laws_concept.png', dpi=140)
print('roots at R=10%%, eps=1 pm: Delta_lo=%.1f pm, Delta_hi=%.1f pm'
      % (dlo, dhi))
print('saved figs/fig_s26_laws_concept.png')
