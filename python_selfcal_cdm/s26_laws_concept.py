"""s26_laws_concept.py - physical meaning of Rule A and Rule B.

Row (a) starts with a sketch of the light path: the launched light crosses
the upstream grating j twice, once toward grating k and once back, so the
line of k reaches the detector multiplied by the two-pass transmission
(1-R_j)^2. The three panels then show that notch in its three regimes:
co-tuned, at the worst detuning sigma*sqrt(3/2), and far away. Each panel
has the undistorted line, the two-pass transmission, the line at the
detector, and the Gaussian fitted to it, so the bias is the visible
distance between two centres, and the worst panel marks the direction:
away from the neighbour, sign opposite to the detuning.

Panel (b) turns the same pairwise bias into a placement rule over the
complete sensor operating range. Panel (c) paints that rule on the address
plane of one grating, in the style of the CDM-WDM addressing figure: the
danger zone is a spectral band, it covers every upstream grating that falls
inside it whatever its distance along the fiber, and nothing downstream of k
shadows k. Panel (d) shows why the mean code-leakage bias behaves like a
wavelength-axis stretch and why two stabilized references remove its offset
and scale.
"""
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
import numpy as np
from scipy.optimize import brentq, curve_fit

import figstyle as FS
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

warnings.filterwarnings('ignore')
FS.apply(base=7.1)

SIG = 106.2
R = 0.10
EPS = 1.0
C_A = (4.0 / 3.0) * np.sqrt(2.0 / 3.0)
DSTAR = SIG * np.sqrt(1.5)
DETS = [(0.0, 'co-tuned, $\\Delta\\lambda_{jk}=0$'),
        (DSTAR, 'worst, $\\Delta\\lambda_{jk}=\\sigma\\sqrt{3/2}$'),
        (400.0, 'far, $\\Delta\\lambda_{jk}=400$ pm')]


def panel_title(ax, letter, text):
    ax.set_title('(%s) %s' % (letter, text), loc='left', pad=4,
                 fontsize=7.6)


def gaussian(x, amp, centre, width, baseline):
    return baseline + amp * np.exp(-0.5 * ((x - centre) / width) ** 2)


fig = plt.figure(figsize=(7.1, 4.15))
gs0 = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.12], hspace=0.55,
                       left=0.052, right=0.995, bottom=0.095, top=0.955)
top = gs0[0].subgridspec(1, 1)
bot = gs0[1].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.28)
pa = fig.add_subplot(top[0, 0])
b = fig.add_subplot(bot[0, 0])
c = fig.add_subplot(bot[0, 1])

# ---------------------------------------------------------------------------
# (a) the three regimes merged on one axis, the light-path sketch as an
# inset in the empty left part (x extended to -730 pm to make that room)
# ---------------------------------------------------------------------------
panel_title(pa, 'a', 'Rule A on the spectrum')
nu = np.linspace(-330, 560, 1800)
wanted = np.exp(-0.5 * (nu / SIG) ** 2)
CASES = ((0.0, FS.BLUE, 'co-tuned'), (DSTAR, FS.VERM, 'worst'),
         (400.0, FS.GREEN, 'far'))
ANN = ((-262, 0.86), (235, 0.66), (455, 0.70))
shifts = []
for (det, colr, name), (tx, ty) in zip(CASES, ANN):
    two_pass = (1.0 - R * np.exp(-0.5 * ((nu - det) / SIG) ** 2)) ** 2
    received = wanted * two_pass
    (amp, mu, sig_, base), _ = curve_fit(
        gaussian, nu, received, p0=[0.9, 0.0, SIG, 0.0], maxfev=20000)
    shifts.append(mu)
    pa.plot(nu, two_pass, color=colr, lw=0.85, alpha=0.85, zorder=2)
    pa.fill_between(nu, received, wanted, where=wanted >= received,
                    color=colr, alpha=0.16, lw=0, zorder=1)
    pa.plot(nu, received, color=colr, lw=1.45, zorder=4)
    pa.axvline(mu, color=colr, lw=0.7, ls=(0, (1.5, 2)), zorder=1)
    lab = ('%.1f' % mu).replace('-0.0', '0.0')
    pa.text(tx, ty, name + '\n' + r'$\delta\lambda_{k\leftarrow j}=%s$ pm' % lab,
            fontsize=5.5, color=colr, ha='center', va='center')
pa.plot(nu, wanted, color='0.12', ls=(0, (4, 2.5)), lw=1.3, zorder=6)
pa.set_xlim(-580, 560)
pa.set_ylim(0.0, 1.12)
pa.set_yticks([0, 0.5, 1.0])
pa.set_xticks([-400, -200, 0, 200, 400])
pa.set_xlabel(r'wavelength offset from $\lambda_{B,k}$ [pm]')
pa.set_ylabel('normalized reflectance')
roles = [Line2D([], [], color='0.3', lw=0.85, alpha=0.85),
         Line2D([], [], color='0.12', ls=(0, (4, 2.5)), lw=1.3),
         Line2D([], [], color='0.3', lw=1.45)]
pa.legend(roles, [r'two-pass transmission $(1-R_j)^2$',
                  r'undistorted $R_k(\lambda)$',
                  r'at the detector, $A_k$'],
          fontsize=5.2, loc='lower right', handlelength=1.5, borderaxespad=0.3)

sk = pa.inset_axes([0.015, 0.06, 0.235, 0.50])
sk.set_xlim(0, 10)
sk.set_ylim(0, 5)
sk.axis('off')
sk.plot([0.3, 9.7], [2.6, 2.6], color='black', lw=4.0, solid_capstyle='butt',
        zorder=1)
for x0, colr, lab in ((3.6, FS.BLUE, 'FBG$_j$'), (7.0, FS.ORANGE, 'FBG$_k$')):
    sk.add_patch(Rectangle((x0, 1.7), 1.3, 1.8, facecolor=colr,
                           edgecolor='none', zorder=2))
    sk.text(x0 + 0.65, 2.62, lab, ha='center', va='center', fontsize=5.6,
            color='white', zorder=3)
sk.annotate('', xy=(6.9, 3.05), xytext=(0.9, 3.05),
            arrowprops=dict(arrowstyle='-|>', lw=0.9, color='0.25'))
sk.text(0.95, 3.34, r'$\times(1{-}R_j)$', fontsize=5.2, color=FS.BLUE)
sk.annotate('', xy=(0.9, 2.15), xytext=(6.9, 2.15),
            arrowprops=dict(arrowstyle='-|>', lw=0.9, color='0.25'))
sk.text(0.95, 1.34, r'$\times(1{-}R_j)$', fontsize=5.2, color=FS.BLUE)
sk.text(7.15, 1.34, r'$\times R_k$', fontsize=5.2, color='#B87F00')
sk.text(5.0, 0.34, r'$A_k=R_k\,(1-R_j)^2$', ha='center', fontsize=5.6,
        color='0.15')
print('Fig. 3(a): fitted shifts %.2f, %.2f, %.2f pm at detunings 0, %.0f, '
      '400 pm' % (shifts[0], shifts[1], shifts[2], DSTAR))

# ---------------------------------------------------------------------------
# (b) Pairwise bias as a placement rule over the full sensor range
# ---------------------------------------------------------------------------
D = np.linspace(0, 700, 1000)
bias = C_A * R * D * np.exp(-D ** 2 / (3.0 * SIG ** 2))
root_fn = lambda x: C_A * R * x * np.exp(-x ** 2 / (3.0 * SIG ** 2)) - EPS
dlo = brentq(root_fn, 1.0, DSTAR)
dhi = brentq(root_fn, DSTAR, 700.0)
half_fn = (lambda x: C_A * R * x * np.exp(-x ** 2 / (3.0 * SIG ** 2))
           - 0.5 * bias.max())
hlo = brentq(half_fn, 1.0, DSTAR)
hhi = brentq(half_fn, DSTAR, 700.0)

b.axvspan(dlo, dhi, color=FS.VERM, alpha=0.13, lw=0)
b.plot(D, bias, color=FS.VERM, lw=1.35, label='Gaussian fit')
bias_cen = (R / np.sqrt(2.0)) * D * np.exp(-D ** 2 / (4.0 * SIG ** 2))
b.plot(D, bias_cen, color=FS.BLUE, lw=1.1, ls=(0, (5, 2)), label='centroid')
b.axhline(EPS, color='0.35', lw=0.75, ls=(0, (4, 2)),
          label='tolerance $\\epsilon=1$ pm')
b.legend(loc='upper right', bbox_to_anchor=(0.985, 0.295), fontsize=5.4,
         handlelength=1.6, labelspacing=0.2, borderaxespad=0.0)
b.axvline(dlo, color=FS.VERM, lw=0.55, ls=(0, (2, 2)))
b.axvline(dhi, color=FS.VERM, lw=0.55, ls=(0, (2, 2)))
b.plot(DSTAR, bias.max(), 'o', color=FS.VERM, ms=3.5)
b.text(150, 9.48,
       'maximum $0.81R\\sigma$\nat $\\sigma\\sqrt{3/2}$',
       fontsize=5.4, color=FS.VERM, va='top', ha='left')

forb = (dlo, dhi)
safe = (dhi + 18.0, 690.0)
b.plot(forb, [-0.62, -0.62], color=FS.VERM, lw=8.0,
       solid_capstyle='butt', clip_on=False)
b.text(np.mean(forb), -0.62, 'forbidden', ha='center', va='center',
       fontsize=5.8, color='white', clip_on=False)
b.plot(safe, [-0.62, -0.62], color=FS.GREEN, lw=8.0,
       solid_capstyle='butt', clip_on=False)
b.text(np.mean(safe), -0.62, 'allowed', ha='center', va='center',
       fontsize=5.8, color='white', clip_on=False)

b.set_xlim(0, 700)
b.set_ylim(-1.15, 9.55)
b.set_yticks([0, 4, 8])
b.set_xlabel(r'pair detuning $|\Delta\lambda_{jk}|$ [pm]')
b.set_ylabel(r'pairwise bias $|\delta\lambda_{k\leftarrow j}|$ [pm]')
panel_title(b, 'b', 'Rule A as a placement criterion')


# ---------------------------------------------------------------------------
# (c) who is priced by Rule A: the address plane of one grating
# ---------------------------------------------------------------------------
XK = 80.0
YLIM = 430.0
m = b.inset_axes([0.50, 0.42, 0.48, 0.52])
m.set_xlim(0, 127)
m.set_ylim(-YLIM, YLIM)
m.set_xticks([0, 80])
m.set_yticks([-400, 0, 400])
m.tick_params(labelsize=5.0, length=2.6, width=1.0, pad=1.5)
m.set_xlabel(r'delay bin $\tau$', fontsize=5.2, labelpad=1.0)
m.set_ylabel(r'$\Delta\lambda_{jk}$ [pm]', fontsize=5.2, labelpad=0.5)

m.axvline(XK, color='0.45', lw=0.6, ls=(0, (2, 2)))
for y0 in (dlo, -dhi):
    m.add_patch(Rectangle((0.0, y0), XK, dhi - dlo,
                          color=FS.VERM, alpha=0.13, lw=0))
m.plot(XK, 0.0, 'v', ms=4.0, color=FS.ORANGE)
m.text(XK + 4, -150.0, '$k$, read', ha='left', fontsize=4.6, color='#B87F00')

up_hit = [(10.0, 150.0), (33.0, -230.0), (58.0, 305.0)]
up_safe = [(20.0, 4.0), (44.0, -380.0)]
down = [(97.0, 140.0), (108.0, -260.0), (120.0, 30.0)]
m.plot([x for x, _ in up_hit], [y for _, y in up_hit], 'v', ms=3.4,
       color=FS.VERM, ls='none')
m.plot([x for x, _ in up_safe + down], [y for _, y in up_safe + down],
       'v', ms=3.4, color=FS.GREEN, ls='none')
m.text(38.0, 358.0, 'shadows $k$', ha='center', fontsize=4.6, color=FS.VERM,
       bbox=dict(fc='white', ec='none', pad=0.4, alpha=0.85))
m.text(100.0, -365.0, 'downstream: safe', ha='center', va='center',
       fontsize=4.6, color=FS.GREEN)

# ---------------------------------------------------------------------------
# (d) Rule B is an axis stretch, removed by two reference anchors
# ---------------------------------------------------------------------------
for sp in ('top', 'right', 'left'):
    c.spines[sp].set_visible(False)
c.set_yticks([])
c.set_xticks([-180, -90, 0, 90, 180])
c.tick_params(axis='x', labelsize=6.0, length=2.2, width=0.6)
c.set_xlabel(r'position in the band $\nu_k=\lambda_{B,k}-\lambda_0$ [pm]')
c.set_xlim(-300, 275)
c.set_ylim(-0.32, 3.08)
panel_title(c, 'c', 'Rule B: an axis stretch, two references remove it')

W = 180.0
K, N = 32, 127
stretch = 1.0 + (K - 1) / N
true = np.linspace(-W, W, 7)
read = true * stretch
rows = [(2.30, true, 'true', '0.35'),
        (0.86, read, 'read, Rule B', FS.VERM)]

for y0, values, label, col in rows:
    c.plot([-235, 235], [y0, y0], color='0.78', lw=0.75)
    for q, x0 in enumerate(values):
        is_ref = q in (0, len(values) - 1)
        tick_col = FS.BLUE if is_ref else col
        c.plot([x0, x0], [y0 - 0.16, y0 + 0.16], color=tick_col,
               lw=2.1 if is_ref else 1.15)
    c.text(-262, y0, label, ha='right', va='center', fontsize=6.2,
           color=col)

for x0, xr in zip(true, read):
    if abs(x0) < 1.0:
        continue
    c.add_patch(FancyArrowPatch((x0, 2.10), (xr, 1.06),
                                arrowstyle='-|>', mutation_scale=5,
                                color=FS.VERM, lw=0.6, alpha=0.75))
c.text(-W, 2.72, 'ref.', ha='center', fontsize=5.8, color=FS.BLUE)
c.text(W, 2.72, 'ref.', ha='center', fontsize=5.8, color=FS.BLUE)
c.text(0.0, 1.56, 'every grating pulled outward by $1+(K-1)/N$',
       ha='center', fontsize=5.7, color=FS.VERM,
       bbox=dict(fc='white', ec='none', pad=0.6, alpha=0.9))
c.text(0.0, 0.26,
       'two anchored references (blue) measure the stretch and take it out',
       ha='center', fontsize=5.4, color=FS.BLUE)

fig.savefig('figs/fig_s26_laws_concept.pdf', bbox_inches='tight',
            pad_inches=0.025)
fig.savefig('figs/fig_s26_laws_concept.png', dpi=300, bbox_inches='tight',
            pad_inches=0.025)
plt.close(fig)

print('Rule A roots at 1 pm: Delta_lo=%.1f pm, Delta_hi=%.1f pm' % (dlo, dhi))
print('Rule A at least half its maximum between %.0f and %.0f pm '
      '(%.2f and %.2f sigma)' % (hlo, hhi, hlo / SIG, hhi / SIG))
print('saved figs/fig_s26_laws_concept.pdf and .png')
