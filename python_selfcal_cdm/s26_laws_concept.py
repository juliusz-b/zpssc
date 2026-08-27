"""s26_laws_concept.py - concept panels for the closed-form laws, v2.

One-glance versions:
  (a) Law A: the neighbour's notch BITES a piece out of the flank of the
      line being read; the bite is shaded, the reading moves by delta.
      No transmission curve, no legend - direct labels only.
  (b) The placement map: |delta(Delta)| with the forbidden band; below the
      axis, two operating windows of a +/-10 K pair with a cross (co-tuned,
      enters the band) and a check (far-detuned, clears it).
  (c) Law B: three rulers. True positions; as read (stretched by
      1+(K-1)/N, true scale); after the two-reference fit (pinned back).
      Connector fans make the stretch visible instantly.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import figstyle as FS

FS.apply()

SIG = 106.2          # pm
C_A = (4.0 / 3.0) * np.sqrt(2.0 / 3.0)

fig, ax = plt.subplots(1, 3, figsize=(7.1, 2.25))

# ------------------------- (a) Law A: the bite ------------------------------
a = ax[0]
nu = np.linspace(-420, 420, 900)
Rx = 0.35                                  # exaggerated for visibility
DET = 130.0
g = np.exp(-0.5 * (nu / SIG) ** 2)
T = (1.0 - Rx * np.exp(-0.5 * ((nu - DET) / SIG) ** 2)) ** 2
seen = g * T
a.fill_between(nu / 1000, g, color='0.88', lw=0)
a.plot(nu / 1000, g, color='0.55', lw=0.9)
a.fill_between(nu / 1000, seen, g, color=FS.BLUE, alpha=0.30, lw=0)
a.plot(nu / 1000, seen, color=FS.VERM, lw=1.4)
mu_seen = nu[np.argmax(seen)] / 1000
a.annotate('', xy=(mu_seen, 1.055), xytext=(0, 1.055),
           arrowprops=dict(arrowstyle='-|>', lw=1.3, color=FS.VERM))
a.plot([0, 0], [0, 1.0], color='0.45', lw=0.7, ls=(0, (2, 2)))
a.plot([mu_seen, mu_seen], [0, seen.max()], color=FS.VERM, lw=0.7,
       ls=(0, (2, 2)))
a.text(-0.028, 1.10, r'reading moves by $\delta$', fontsize=7.0,
       color=FS.VERM, ha='center')
a.text(0.21, 0.66, 'bite by the\nneighbour\n' r'detuned by $\Delta$',
       fontsize=6.6, color=FS.BLUE, ha='center')
a.annotate('', xy=(0.115, 0.72), xytext=(0.175, 0.72),
           arrowprops=dict(arrowstyle='-|>', lw=0.6, color=FS.BLUE))
a.text(-0.27, 0.80, 'true\nline', fontsize=6.6, color='0.4', ha='center')
a.text(-0.115, 0.38, 'line\nas read', fontsize=6.6, color=FS.VERM,
       ha='center')
a.set_xlim(-0.42, 0.42); a.set_ylim(0, 1.22)
a.set_yticks([])
a.set_xlabel(r'$\lambda$ offset [nm]')
a.set_title('(a) Law A: the notch pulls the reading', fontsize=9)
for sp in ('left',):
    a.spines[sp].set_visible(False)

# ------------------------- (b) the placement map ----------------------------
b = ax[1]
R = 0.10
eps = 1.0
D = np.linspace(0, 700, 800)
bias = C_A * R * D * np.exp(-D ** 2 / (3 * SIG ** 2))
from scipy.optimize import brentq
f = lambda x: C_A * R * x * np.exp(-x ** 2 / (3 * SIG ** 2)) - eps
dlo = brentq(f, 1, SIG * np.sqrt(1.5))
dhi = brentq(f, SIG * np.sqrt(1.5), 700)
b.axvspan(dlo, dhi, color=FS.VERM, alpha=0.12, lw=0)
b.plot(D, bias, color=FS.VERM, lw=1.3)
b.axhline(eps, color='0.35', lw=0.7, ls=(0, (4, 2)))
b.text(690, 1.28, r'tolerance $\epsilon = 1$ pm', fontsize=6.2, ha='right',
       color='0.35')
b.text((dlo + dhi) / 2, 8.55, 'forbidden', fontsize=7.5, ha='center',
       color=FS.VERM)
for dv, lab, ha in [(dlo + 14, r'$\Delta_\mathrm{lo}$', 'left'),
                    (dhi, r'$\Delta_\mathrm{hi}$', 'center')]:
    b.text(dv, 9.10, lab, fontsize=7, ha=ha, color=FS.VERM)
for dv in (dlo, dhi):
    b.axvline(dv, color=FS.VERM, lw=0.6, ls=(0, (2, 2)))
# operating windows strip under the curve
S = 100.0
b.plot([0, 2 * S], [-1.75, -1.75], color='0.2', lw=4,
       solid_capstyle='butt', clip_on=False)
b.text(2 * S + 25, -1.75, r'$\times$', fontsize=10, color=FS.VERM,
       va='center')
b.text(2 * S + 75, -1.75, 'co-tuned pair: biased at some $T$',
       fontsize=6.2, va='center', color='0.2')
b.plot([550 - 2 * S, 550 + 2 * S], [-3.25, -3.25], color=FS.GREEN, lw=4,
       solid_capstyle='butt', clip_on=False)
b.text(550 + 2 * S + 25, -3.25, r'$\checkmark$', fontsize=10,
       color=FS.GREEN, va='center')
b.text(550 - 2 * S - 25, -3.25, 'far-detuned: safe at every $T$',
       fontsize=6.2, va='center', ha='right', color=FS.GREEN)
b.set_xlim(0, 700); b.set_ylim(0, 9.8)
b.set_yticks([0, 4, 8])
b.set_xlabel(r'pair detuning $\Delta$ [pm]', labelpad=26)
b.set_ylabel(r'pairwise bias $|\delta|$ [pm]')
b.set_title(r'(b) windows of a $\pm10$ K pair', fontsize=9)

# ------------------------- (c) Law B: three rulers --------------------------
c = ax[2]
W = 200.0
pos = np.linspace(-W, W, 9)
K, N = 32, 127
st = 1.0 + (K - 1) / N                     # true stretch, 1.244
rows = [('true positions', pos, '0.45', 2.0),
        ('as read: stretched by $1+(K{-}1)/N$', pos * st, FS.VERM, 1.0),
        ('after the 2-reference fit', pos, FS.BLUE, 0.0)]
REF = [-150, 150]
for lab, xs, col, yy in rows:
    c.plot([-W * st - 20, W * st + 20], [yy, yy], color='0.85', lw=0.8)
    for x0, xt in zip(pos, xs):
        isref = any(abs(x0 - r0) < 1 for r0 in REF)
        cc = FS.BLUE if isref else col
        c.plot([xt, xt], [yy - 0.13, yy + 0.13], color=cc,
               lw=2.2 if isref else 1.1)
    c.text(-W * st - 35, yy, lab, fontsize=6.4, ha='right', va='center',
           color=col if col != '0.45' else '0.3')
# connector fan true -> stretched
for x0 in pos:
    c.plot([x0, x0 * st], [2.0 - 0.16, 1.0 + 0.16], color='0.75', lw=0.5)
# connector fan stretched -> corrected (pinned at refs)
for x0 in pos:
    c.plot([x0 * st, x0], [1.0 - 0.16, 0.0 + 0.16], color=FS.BLUE,
           lw=0.5, alpha=0.55)
c.text(0, 2.62, r'$9.4 \rightarrow 2.3$ pm at $K=32$, $N=127$',
       fontsize=6.6, ha='center', color='0.15')
c.set_xlim(-W * st - 150, W * st + 40)
c.set_ylim(-0.45, 2.95)
c.set_xticks([]); c.set_yticks([])
for sp in c.spines.values():
    sp.set_visible(False)
c.set_xlabel(r'position in the band $\nu_k$')
c.set_title('(c) Law B: a stretch, pinned by references', fontsize=9)

plt.tight_layout()
fig.savefig('figs/fig_s26_laws_concept.png', dpi=140,
            bbox_inches='tight', pad_inches=0.03)
print('roots: Delta_lo=%.1f, Delta_hi=%.1f pm' % (dlo, dhi))
print('saved figs/fig_s26_laws_concept.png')
