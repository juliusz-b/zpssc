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
                       left=0.052, right=0.995, bottom=0.095, top=0.92)
top = gs0[0].subgridspec(1, 4, width_ratios=[0.68, 1.0, 1.0, 1.0],
                         wspace=0.24)
bot = gs0[1].subgridspec(1, 3, width_ratios=[0.98, 0.92, 1.10], wspace=0.36)
sk = fig.add_subplot(top[0, 0])
axa = [fig.add_subplot(top[0, i]) for i in (1, 2, 3)]
b = fig.add_subplot(bot[0, 0])
m = fig.add_subplot(bot[0, 1])
c = fig.add_subplot(bot[0, 2])

# ---------------------------------------------------------------------------
# (a0) the light path: two passes through j, one reflection at k
# ---------------------------------------------------------------------------
sk.set_xlim(0, 10)
sk.set_ylim(0, 10)
sk.axis('off')
panel_title(sk, 'a', 'Rule A on the spectrum')

sk.plot([0.4, 9.6], [5.3, 5.3], color='0.55', lw=1.6)
for dx in (-0.22, 0.0, 0.22):
    sk.plot([3.3 + dx, 3.3 + dx], [4.5, 6.1], color=FS.BLUE, lw=1.3)
    sk.plot([7.3 + dx, 7.3 + dx], [4.5, 6.1], color=FS.ORANGE, lw=1.3)
sk.text(3.3, 6.6, '$j$, upstream', ha='center', fontsize=5.8, color=FS.BLUE)
sk.text(7.62, 6.6, '$k$, read', ha='left', fontsize=5.8, color=FS.ORANGE)

# launch lane, reflection hook at k, return lane
sk.annotate('', xy=(7.15, 7.8), xytext=(0.5, 7.8),
            arrowprops=dict(arrowstyle='-|>', color='0.35', lw=0.9,
                            mutation_scale=7))
sk.plot([7.3, 7.3], [7.8, 3.0], color='0.35', lw=0.8)
sk.annotate('', xy=(0.5, 3.0), xytext=(7.3, 3.0),
            arrowprops=dict(arrowstyle='-|>', color='0.35', lw=0.9,
                            mutation_scale=7))
sk.text(3.3, 8.3, r'$\times(1{-}R_j)$', ha='center', fontsize=5.6,
        color=FS.BLUE)
sk.text(3.3, 2.0, r'$\times(1{-}R_j)$', ha='center', fontsize=5.6,
        color=FS.BLUE)
sk.text(7.75, 4.0, r'$\times R_k$', ha='left', fontsize=5.6,
        color=FS.ORANGE)

sk.text(5.0, 0.6, r'$A_k=R_k\,(1-R_j)^2$', ha='center', fontsize=6.6,
        color='0.15')

# ---------------------------------------------------------------------------
# (a1-a3) the three regimes, on the spectrum itself
# ---------------------------------------------------------------------------
nu = np.linspace(-5.0 * SIG, 5.0 * SIG, 1800)
wanted = np.exp(-0.5 * (nu / SIG) ** 2)
shifts = []
for a, (det, label) in zip(axa, DETS):
    upstream = R * np.exp(-0.5 * ((nu - det) / SIG) ** 2)
    two_pass = (1.0 - upstream) ** 2
    received = wanted * two_pass
    (amp, mu, sig, base), _ = curve_fit(
        gaussian, nu, received, p0=[0.9, 0.0, SIG, 0.0], maxfev=20000)
    fit = gaussian(nu, amp, mu, sig, base)
    shifts.append(mu)

    a.plot(nu, two_pass, color=FS.BLUE, lw=0.8,
           label='two-pass transmission $(1-R_j)^2$')
    a.plot(nu, wanted, color='0.50', lw=1.0, ls=(0, (3, 2)),
           label='undistorted $R_k(\\lambda)$')
    a.fill_between(nu, received, wanted, where=wanted >= received,
                   color=FS.ORANGE, alpha=0.28, lw=0)
    a.plot(nu, received, color=FS.ORANGE, lw=1.25, label='at the detector, $A_k$')
    a.plot(nu, fit, color=FS.VERM, lw=0.9, ls=(0, (3, 1.7)),
           label='Gaussian fit $\\to\\hat\\lambda_{B,k}$')
    a.axvline(0.0, color='0.45', lw=0.65, ls=(0, (2, 2)))
    a.axvline(mu, color=FS.VERM, lw=0.75, ls=(0, (2, 2)))
    lab = ('%.1f' % mu).replace('-0.0', '0.0')
    if abs(mu) > 1.0:
        FS.dim_gap(a, 0.0, mu, 1.13,
                   r'$\delta\lambda_{k\leftarrow j}=%s$ pm' % lab,
                   color=FS.VERM, tail=55.0, side='right', fontsize=6.0)
        a.annotate('away from $j$', xy=(mu - 14, 0.90),
                   xytext=(-320, 0.66), fontsize=5.6, color=FS.VERM,
                   va='center',
                   arrowprops=dict(arrowstyle='-|>', color=FS.VERM, lw=0.6,
                                   mutation_scale=6))
    else:
        a.text(8.0, 1.13, r'$\delta\lambda_{k\leftarrow j}=%s$ pm' % lab,
               fontsize=6.0, color=FS.VERM, va='center')
    if det > 1.0:
        a.text(det, 0.71, '$j$', ha='center', fontsize=6.2, color=FS.BLUE)
    a.text(0.975, 0.05, label, transform=a.transAxes, ha='right',
           va='bottom', fontsize=5.9, color='0.25')
    a.set_xlim(-330, 500)
    a.set_ylim(0.0, 1.30)
    a.set_yticks([0, 0.5, 1.0])
    a.set_xticks([-200, 0, 200, 400])
    a.set_xlabel(r'wavelength offset from $\lambda_{B,k}$ [pm]')
axa[0].set_ylabel('normalized reflectance')
for a in axa[1:]:
    a.set_yticklabels([])
handles, labels = axa[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.995, 1.005),
           ncol=4, fontsize=5.6, frameon=False, handlelength=1.4,
           columnspacing=0.9, borderaxespad=0.1)
print('Fig. 2(a): fitted shifts %.2f, %.2f, %.2f pm at detunings 0, %.0f, '
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
b.plot(D, bias, color=FS.VERM, lw=1.35, label='Rule A, $|\\delta\\lambda_{k\\leftarrow j}|$')
b.axhline(EPS, color='0.35', lw=0.75, ls=(0, (4, 2)),
          label='tolerance $\\epsilon=1$ pm')
b.axvline(dlo, color=FS.VERM, lw=0.55, ls=(0, (2, 2)))
b.axvline(dhi, color=FS.VERM, lw=0.55, ls=(0, (2, 2)))
b.plot(DSTAR, bias.max(), 'o', color=FS.VERM, ms=3.5)
b.text(DSTAR + 14, bias.max() + 0.05,
       'maximum $0.81R\\sigma$\nat $\\sigma\\sqrt{3/2}$',
       fontsize=5.6, color=FS.VERM, va='center', ha='left')

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
b.legend(loc='upper right', fontsize=5.7, frameon=False, handlelength=1.5,
         labelspacing=0.18, borderaxespad=0.25)

# ---------------------------------------------------------------------------
# (c) who is priced by Rule A: the address plane of one grating
# ---------------------------------------------------------------------------
XK = 80.0                       # delay bin of the grating being read
YLIM = 430.0
m.set_xlim(0, 127)
m.set_ylim(-YLIM, YLIM)
m.set_xticks([0, 40, 80, 120])
m.set_yticks([-400, -200, 0, 200, 400])
m.set_xlabel(r'delay bin $\tau$ (position along fiber)')
m.set_ylabel(r'detuning $\Delta\lambda_{jk}$ [pm]')
panel_title(m, 'c', 'the pairs Rule A prices')

m.axvline(XK, color='0.45', lw=0.7, ls=(0, (2, 2)))
for y0 in (dlo, -dhi):
    m.add_patch(Rectangle((0.0, y0), XK, dhi - dlo,
                          color=FS.VERM, alpha=0.13, lw=0))
    m.plot([0.0, XK], [y0, y0], color=FS.VERM, lw=0.55, ls=(0, (2, 2)))
    m.plot([0.0, XK], [y0 + dhi - dlo, y0 + dhi - dlo],
           color=FS.VERM, lw=0.55, ls=(0, (2, 2)))

m.plot(XK, 0.0, 'v', ms=6.0, color=FS.ORANGE, zorder=5)
m.text(XK, -95.0, '$k$, read', ha='center', fontsize=5.8, color=FS.ORANGE)

up_hit = [(10.0, 150.0), (33.0, -230.0), (58.0, 305.0)]
up_safe = [(20.0, 0.0), (44.0, -395.0)]
down = [(97.0, 140.0), (108.0, -260.0), (120.0, 30.0)]
for x0, y0 in up_hit:
    m.plot(x0, y0, 'v', ms=5.0, color=FS.VERM)
for x0, y0 in up_safe + down:
    m.plot(x0, y0, 'v', ms=5.0, color=FS.GREEN)

m.text(40.0, YLIM * 0.86, 'upstream of $k$', ha='center', fontsize=5.8,
       color='0.30')
m.text(104.0, YLIM * 0.86, 'downstream,\nnever shadows $k$', ha='center',
       fontsize=5.8, color=FS.GREEN)
m.text(34.0, 232.0, 'biased pairs,\nany distance', ha='center', fontsize=5.6,
       color=FS.VERM)
m.text(20.0, -78.0, 'co-tuned, safe', ha='center', fontsize=5.4,
       color=FS.GREEN, bbox=dict(fc='white', ec='none', pad=0.5, alpha=0.85))
m.text(49.0, -400.0, 'far, safe', ha='left', va='center', fontsize=5.4,
       color=FS.GREEN)

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
panel_title(c, 'd', 'Rule B: an axis stretch, two references remove it')

W = 180.0
K, N = 32, 127
stretch = 1.0 + (K - 1) / N
true = np.linspace(-W, W, 7)
read = true * stretch
rows = [(2.38, true, 'true', '0.35'),
        (1.43, read, 'read, Rule B', FS.VERM),
        (0.48, true, 'calibrated', FS.BLUE)]

for y0, values, label, col in rows:
    c.plot([-205, 205], [y0, y0], color='0.78', lw=0.75)
    for q, x0 in enumerate(values):
        is_ref = q in (0, len(values) - 1)
        tick_col = FS.BLUE if is_ref else col
        c.plot([x0, x0], [y0 - 0.13, y0 + 0.13], color=tick_col,
               lw=2.1 if is_ref else 1.15)
    c.text(-265, y0, label, ha='right', va='center', fontsize=6.2,
           color=col)

for x0, xr in zip(true, read):
    c.add_patch(FancyArrowPatch((x0, 2.22), (xr, 1.59),
                                arrowstyle='-', color=FS.VERM,
                                lw=0.55, alpha=0.55))
    c.add_patch(FancyArrowPatch((xr, 1.27), (x0, 0.64),
                                arrowstyle='-', color=FS.BLUE,
                                lw=0.55, alpha=0.55))

c.text(-W, 2.70, 'ref.', ha='center', fontsize=5.8, color=FS.BLUE)
c.text(W, 2.70, 'ref.', ha='center', fontsize=5.8, color=FS.BLUE)
c.text(0.0, 1.90, 'every grating pulled outward by $1+(K-1)/N$',
       ha='center', fontsize=5.6, color=FS.VERM,
       bbox=dict(fc='white', ec='none', pad=0.6, alpha=0.85))

fig.savefig('figs/fig_s26_laws_concept.pdf', bbox_inches='tight',
            pad_inches=0.025)
fig.savefig('figs/fig_s26_laws_concept.png', dpi=300, bbox_inches='tight',
            pad_inches=0.025)
plt.close(fig)

print('Rule A roots at 1 pm: Delta_lo=%.1f pm, Delta_hi=%.1f pm' % (dlo, dhi))
print('Rule A at least half its maximum between %.0f and %.0f pm '
      '(%.2f and %.2f sigma)' % (hlo, hhi, hlo / SIG, hhi / SIG))
print('saved figs/fig_s26_laws_concept.pdf and .png')
