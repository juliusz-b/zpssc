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


fig = plt.figure(figsize=(7.1, 2.55))
gs0 = fig.add_gridspec(1, 2, width_ratios=[1.55, 1.0], wspace=0.30,
                       left=0.052, right=0.995, bottom=0.17, top=0.90)
pa = fig.add_subplot(gs0[0, 0])
b = fig.add_subplot(gs0[0, 1])

# ---------------------------------------------------------------------------
# (a) the three regimes merged on one axis, the light-path sketch as an
# inset in the empty left part (x extended to -730 pm to make that room)
# ---------------------------------------------------------------------------
panel_title(pa, 'a', 'Rule A on the spectrum')
nu = np.linspace(-330, 560, 1800)
wanted = np.exp(-0.5 * (nu / SIG) ** 2)
CASES = ((0.0, FS.BLUE, 'co-tuned'), (DSTAR, FS.VERM, 'worst'),
         (400.0, FS.GREEN, 'far'))
ANN = ((-250, 0.62), (235, 0.66), (455, 0.70))
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
          fontsize=5.0, loc='upper left', bbox_to_anchor=(0.005, 0.985),
          handlelength=1.5, borderaxespad=0.0)

sk = pa.inset_axes([0.015, 0.06, 0.27, 0.52])
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
b.plot(D, bias, color=FS.VERM, lw=1.35, label='Rule A, $\\delta\\lambda_{k\\leftarrow j}$')
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
ins = b.inset_axes([0.40, 0.44, 0.58, 0.50])
ins.set_xlim(0, 10)
ins.set_ylim(0, 5)
ins.axis('off')
ins.plot([0.2, 9.8], [2.5, 2.5], color='black', lw=3.0, solid_capstyle='butt',
         zorder=1)
EX = [(1.0, FS.VERM, '$j_1$', '+150 pm', 'biases $k$'),
      (3.0, FS.GREEN, '$j_2$', '0 pm', 'safe'),
      (5.0, FS.GREEN, '$j_3$', '+400 pm', 'safe'),
      (7.0, FS.ORANGE, '$k$', 'read', ''),
      (9.0, FS.GREEN, '$d$', 'behind $k$', 'safe')]
for x0, colr, lab, det, verdict in EX:
    ins.add_patch(Rectangle((x0 - 0.42, 1.75), 0.84, 1.5, facecolor=colr,
                            edgecolor='none', zorder=2))
    ins.text(x0, 2.5, lab, ha='center', va='center', fontsize=5.2,
             color='white', zorder=3)
    ins.text(x0, 1.35, det, ha='center', va='top', fontsize=4.6,
             color=('#B87F00' if colr == FS.ORANGE else colr))
    if verdict:
        ins.text(x0, 3.55, verdict, ha='center', va='bottom', fontsize=4.6,
                 color=colr)
ins.annotate('', xy=(9.6, 0.5), xytext=(0.4, 0.5),
             arrowprops=dict(arrowstyle='-|>', lw=0.7, color='0.4',
                             mutation_scale=6))
ins.text(5.0, 0.1, 'from the laser', ha='center', va='top', fontsize=4.6,
         color='0.4')
ins.text(5.0, 4.9, r'example: detuning of each grating from $k$',
         ha='center', va='top', fontsize=4.8, color='0.25')

fig.savefig('figs/fig_s26_laws_concept.pdf', bbox_inches='tight',
            pad_inches=0.025)
fig.savefig('figs/fig_s26_laws_concept.png', dpi=300, bbox_inches='tight',
            pad_inches=0.025)
plt.close(fig)

print('Rule A roots at 1 pm: Delta_lo=%.1f pm, Delta_hi=%.1f pm' % (dlo, dhi))
print('Rule A at least half its maximum between %.0f and %.0f pm '
      '(%.2f and %.2f sigma)' % (hlo, hhi, hlo / SIG, hhi / SIG))
print('saved figs/fig_s26_laws_concept.pdf and .png')
