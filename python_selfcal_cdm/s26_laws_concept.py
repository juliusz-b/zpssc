"""s26_laws_concept.py - physical meaning of Law A and Law B.

Panel (a) follows the downstream spectrum through the two-pass transmission
notch of one upstream grating. Panel (b) turns the same pairwise bias into a
placement rule over the complete sensor operating range. Panel (c) shows why
the mean code-leakage bias behaves like a wavelength-axis stretch and why two
stabilized references remove its offset and scale.
"""
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import numpy as np
from scipy.optimize import brentq, curve_fit

import figstyle as FS

warnings.filterwarnings('ignore')
FS.apply(base=7.1)

SIG = 106.2
R = 0.10
DET = 130.0
EPS = 1.0
C_A = (4.0 / 3.0) * np.sqrt(2.0 / 3.0)


def panel_title(ax, letter, text):
    ax.set_title(r'$\bf{%s}$  %s' % (letter, text), loc='left', pad=4,
                 fontsize=7.6)


fig, ax = plt.subplots(1, 3, figsize=(7.1, 2.62),
                       gridspec_kw={'width_ratios': [1.10, 1.12, 1.18]})

# ---------------------------------------------------------------------------
# (a) One upstream transmission notch distorts the downstream return
# ---------------------------------------------------------------------------
a = ax[0]
nu = np.linspace(-5.0 * SIG, 5.0 * SIG, 1800)
wanted = np.exp(-0.5 * (nu / SIG) ** 2)
upstream = R * np.exp(-0.5 * ((nu - DET) / SIG) ** 2)
two_pass = (1.0 - upstream) ** 2
received = wanted * two_pass
gaussian = lambda x, amp, centre, width, baseline: \
    baseline + amp * np.exp(-0.5 * ((x - centre) / width) ** 2)
(amp, mu, sig, base), _ = curve_fit(
    gaussian, nu, received, p0=[0.9, 0.0, SIG, 0.0], maxfev=20000)
fit = base + amp * np.exp(-0.5 * ((nu - mu) / sig) ** 2)

a.plot(nu, wanted, color='0.50', lw=1.0, ls=(0, (3, 2)))
a.fill_between(nu, received, wanted, where=wanted >= received,
               color=FS.ORANGE, alpha=0.28, lw=0)
a.plot(nu, received, color=FS.ORANGE, lw=1.25)
a.plot(nu, fit, color=FS.VERM, lw=0.9, ls=(0, (3, 1.7)))
a.axvline(0.0, color='0.45', lw=0.65, ls=(0, (2, 2)))
a.axvline(mu, color=FS.VERM, lw=0.75, ls=(0, (2, 2)))
a.annotate('', xy=(mu, 1.035), xytext=(0.0, 1.035),
           arrowprops=dict(arrowstyle='<->', color=FS.VERM, lw=0.8))
a.text(0.5 * mu, 1.085, r'fitted shift $\delta$', color=FS.VERM,
       fontsize=6.5, ha='center')
a.annotate(r'upstream centre $\Delta=130$ pm', xy=(DET, 0.46),
           xytext=(235, 0.73), fontsize=6.1, color=FS.BLUE, ha='center',
           arrowprops=dict(arrowstyle='-|>', color=FS.BLUE, lw=0.65))
a.text(212, 0.24, 'power removed\nfrom one flank', color=FS.ORANGE,
       fontsize=6.2, ha='center')
a.set_xlim(-320, 340)
a.set_ylim(0.0, 1.18)
a.set_yticks([0, 0.5, 1.0])
a.set_xlabel(r'wavelength offset from $\lambda_{B,k}$ [pm]')
a.set_ylabel('normalized reflectance')
panel_title(a, 'a', 'Law A: a notch moves the fitted centre')
a.plot([-300, -270], [1.12, 1.12], color='0.50', lw=1.0, ls=(0, (3, 2)))
a.text(-262, 1.12, 'wanted line', va='center', fontsize=5.8, color='0.38')
a.plot([-300, -270], [1.05, 1.05], color=FS.ORANGE, lw=1.25)
a.text(-262, 1.05, 'after upstream notch', va='center', fontsize=5.8,
       color=FS.ORANGE)
a.plot([-300, -270], [0.98, 0.98], color=FS.VERM, lw=0.9, ls=(0, (3, 1.7)))
a.text(-262, 0.98, 'Gaussian fit', va='center', fontsize=5.8, color=FS.VERM)

# ---------------------------------------------------------------------------
# (b) Pairwise bias as a placement rule over the full sensor range
# ---------------------------------------------------------------------------
b = ax[1]
D = np.linspace(0, 700, 1000)
bias = C_A * R * D * np.exp(-D ** 2 / (3.0 * SIG ** 2))
root_fn = lambda x: C_A * R * x * np.exp(-x ** 2 / (3.0 * SIG ** 2)) - EPS
dlo = brentq(root_fn, 1.0, SIG * np.sqrt(1.5))
dhi = brentq(root_fn, SIG * np.sqrt(1.5), 700.0)
dstar = SIG * np.sqrt(1.5)

b.axvspan(dlo, dhi, color=FS.VERM, alpha=0.13, lw=0)
b.plot(D, bias, color=FS.VERM, lw=1.35)
b.axhline(EPS, color='0.35', lw=0.75, ls=(0, (4, 2)))
b.axvline(dlo, color=FS.VERM, lw=0.55, ls=(0, (2, 2)))
b.axvline(dhi, color=FS.VERM, lw=0.55, ls=(0, (2, 2)))
b.plot(dstar, bias.max(), 'o', color=FS.VERM, ms=3.5)
b.text(dstar, bias.max() + 0.48, 'worst overlap', ha='center',
       fontsize=6.1, color=FS.VERM)
b.text(0.5 * (dlo + dhi), 6.3, 'forbidden detuning', ha='center',
       fontsize=6.4, color=FS.VERM)
b.text(dlo, 0.25, r'$\Delta_{\rm lo}$', fontsize=6.0, color=FS.VERM,
       ha='right')
b.text(dhi, 0.25, r'$\Delta_{\rm hi}$', fontsize=6.0, color=FS.VERM,
       ha='left')
b.text(690, EPS + 0.25, r'tolerance $\epsilon=1$ pm', ha='right',
       fontsize=6.1, color='0.35')

# Full pair-detuning ranges. The bar, not just its centre, must avoid red.
unsafe = (0.0, 200.0)
safe = (dhi + 18.0, 690.0)
b.plot(unsafe, [-0.82, -0.82], color=FS.VERM, lw=4.0,
       solid_capstyle='butt', clip_on=False)
b.text(np.mean(unsafe), -1.25, 'operating range crosses red', ha='center',
       fontsize=5.9, color=FS.VERM, clip_on=False)
b.plot(safe, [-1.72, -1.72], color=FS.GREEN, lw=4.0,
       solid_capstyle='butt', clip_on=False)
b.text(np.mean(safe), -1.30, 'full range stays safe', ha='center',
       fontsize=5.9, color=FS.GREEN, clip_on=False)

b.set_xlim(0, 700)
b.set_ylim(-2.55, 9.55)
b.set_yticks([0, 4, 8])
b.set_xlabel(r'pair detuning $|\Delta\lambda_{jk}|$ [pm]')
b.set_ylabel(r'pairwise bias $|\delta\lambda_{k\leftarrow j}|$ [pm]')
panel_title(b, 'b', 'Placement: keep the full range outside red')

# ---------------------------------------------------------------------------
# (c) Law B is an axis stretch, removed by two reference anchors
# ---------------------------------------------------------------------------
c = ax[2]
c.axis('off')
c.set_xlim(-335, 275)
c.set_ylim(-0.25, 3.35)
panel_title(c, 'c', 'Law B: two references undo the stretch')

W = 180.0
K, N = 32, 127
stretch = 1.0 + (K - 1) / N
true = np.linspace(-W, W, 7)
read = true * stretch
rows = [(2.55, true, 'true wavelength positions', '0.35'),
        (1.50, read, r'after leakage: $\times[1+(K-1)/N]$', FS.VERM),
        (0.45, true, 'after two-reference calibration', FS.BLUE)]

for y0, values, label, col in rows:
    c.plot([-255, 255], [y0, y0], color='0.78', lw=0.75)
    for q, x0 in enumerate(values):
        is_ref = q in (0, len(values) - 1)
        tick_col = FS.BLUE if is_ref else col
        c.plot([x0, x0], [y0 - 0.13, y0 + 0.13], color=tick_col,
               lw=2.1 if is_ref else 1.15)
    c.text(-270, y0, label, ha='right', va='center', fontsize=6.1,
           color=col)

for x0, xr in zip(true, read):
    c.add_patch(FancyArrowPatch((x0, 2.38), (xr, 1.67),
                                arrowstyle='-', color=FS.VERM,
                                lw=0.55, alpha=0.55))
    c.add_patch(FancyArrowPatch((xr, 1.33), (x0, 0.62),
                                arrowstyle='-', color=FS.BLUE,
                                lw=0.55, alpha=0.55))

c.text(-W, 2.84, 'reference', ha='center', fontsize=5.7, color=FS.BLUE)
c.text(W, 2.84, 'reference', ha='center', fontsize=5.7, color=FS.BLUE)
c.text(0, 3.18, r'mean pull $\propto \nu_k(K-1)/N$ and is independent of $R$',
       ha='center', fontsize=6.2, color='0.25')
c.text(0, -0.02, r'anchors fix offset and scale: $9.4\rightarrow2.3$ pm RMS',
       ha='center', fontsize=6.2, color=FS.BLUE)

fig.subplots_adjust(left=0.075, right=0.995, bottom=0.20, top=0.90,
                    wspace=0.34)
fig.savefig('figs/fig_s26_laws_concept.pdf', bbox_inches='tight',
            pad_inches=0.025)
fig.savefig('figs/fig_s26_laws_concept.png', dpi=300, bbox_inches='tight',
            pad_inches=0.025)
plt.close(fig)

print('Law A roots: Delta_lo=%.1f pm, Delta_hi=%.1f pm' % (dlo, dhi))
print('Panel (a) fitted shift: %.2f pm at Delta=%.0f pm, R=%.0f%%' %
      (mu, DET, 100 * R))
print('saved figs/fig_s26_laws_concept.pdf and .png')
