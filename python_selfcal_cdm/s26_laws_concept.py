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


fig = plt.figure(figsize=(7.1, 2.4))
gs0 = fig.add_gridspec(1, 3, width_ratios=[1.45, 1.0, 1.0], wspace=0.36,
                       left=0.05, right=0.995, bottom=0.18, top=0.90)
pa = fig.add_subplot(gs0[0, 0])
b = fig.add_subplot(gs0[0, 1])
c = fig.add_subplot(gs0[0, 2])

# ---------------------------------------------------------------------------
# (a) the three regimes merged on one axis, the light-path sketch as an
# inset in the empty left part (x extended to -730 pm to make that room)
# ---------------------------------------------------------------------------
panel_title(pa, 'a', 'Rule A on the spectrum')
nu = np.linspace(-330, 560, 1800)
wanted = np.exp(-0.5 * (nu / SIG) ** 2)
CASES = ((0.0, FS.BLUE, 'co-tuned'), (DSTAR, FS.VERM, 'worst'),
         (400.0, FS.GREEN, 'far'))
ANN = ((-215, 0.60), (235, 0.66), (470, 0.50))
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
pa.set_xlim(-330, 560)
pa.set_ylim(0.0, 1.12)
pa.set_yticks([0, 0.5, 1.0])
pa.set_xticks([-200, 0, 200, 400])
pa.set_xlabel(r'wavelength offset from $\lambda_{B,k}$ [pm]')
pa.set_ylabel('normalized reflectance')
roles = [Line2D([], [], color='0.3', lw=0.85, alpha=0.85),
         Line2D([], [], color='0.12', ls=(0, (4, 2.5)), lw=1.3),
         Line2D([], [], color='0.3', lw=1.45)]
pa.legend(roles, [r'$(1-R_j)^2$', r'$R_k(\lambda)$', r'$A_k$'],
          fontsize=4.7, loc='lower right', handlelength=1.2, borderaxespad=0.25,
          handletextpad=0.5)


print('Fig. 3(a): fitted shifts %.2f, %.2f, %.2f pm at detunings 0, %.0f, '
      '400 pm' % (shifts[0], shifts[1], shifts[2], DSTAR))

# ---------------------------------------------------------------------------
# (b) Pairwise bias as a placement rule over the full sensor range
# ---------------------------------------------------------------------------
D = np.linspace(0, 700, 1000)
bias = -C_A * R * D * np.exp(-D ** 2 / (3.0 * SIG ** 2))      # signed, (8)
root_fn = lambda x: C_A * R * x * np.exp(-x ** 2 / (3.0 * SIG ** 2)) - EPS
dlo = brentq(root_fn, 1.0, DSTAR)
dhi = brentq(root_fn, DSTAR, 700.0)
half_fn = (lambda x: C_A * R * x * np.exp(-x ** 2 / (3.0 * SIG ** 2))
           + 0.5 * bias.min())
hlo = brentq(half_fn, 1.0, DSTAR)
hhi = brentq(half_fn, DSTAR, 700.0)

b.axvspan(dlo, dhi, color=FS.VERM, alpha=0.13, lw=0)
b.plot(D, bias, color=FS.VERM, lw=1.35, label='Rule A, (8)')
# the full two-pass model with an unrestricted Gaussian fit, as markers
Dm = np.linspace(25.0, 625.0, 13)
mod = []
for d in Dm:
    rec = wanted * (1.0 - R * np.exp(-0.5 * ((nu - d) / SIG) ** 2)) ** 2
    (amp, mu, sg_, base), _ = curve_fit(
        gaussian, nu, rec, p0=[0.9, 0.0, SIG, 0.0], maxfev=20000)
    mod.append(mu)
b.plot(Dm, mod, 's', ms=3.6, mfc='none', mec=FS.VERM, mew=0.9, ls='none',
       label='full model')
b.axhline(-EPS, color='0.35', lw=0.75, ls=(0, (4, 2)),
          label='tolerance $\\epsilon=1$ pm')
b.axvline(dlo, color=FS.VERM, lw=0.55, ls=(0, (2, 2)))
b.axvline(dhi, color=FS.VERM, lw=0.55, ls=(0, (2, 2)))
b.plot(DSTAR, bias.min(), 'o', color=FS.VERM, ms=3.5)
b.text(150, -9.45, 'min. $-0.81R\\sigma$',
       fontsize=5.0, color=FS.VERM, va='bottom', ha='left')

forb = (dlo, dhi)
safe = (dhi + 18.0, 690.0)
b.plot(forb, [0.62, 0.62], color=FS.VERM, lw=8.0,
       solid_capstyle='butt', clip_on=False)
b.text(np.mean(forb), 0.62, 'forbidden', ha='center', va='center',
       fontsize=5.8, color='white', clip_on=False)
b.plot(safe, [0.62, 0.62], color=FS.GREEN, lw=8.0,
       solid_capstyle='butt', clip_on=False)
b.text(np.mean(safe), 0.62, 'allowed', ha='center', va='center',
       fontsize=5.8, color='white', clip_on=False)

b.set_xlim(0, 700)
b.set_ylim(-9.55, 1.15)
b.set_yticks([-8, -4, 0])
b.set_xlabel(r'pair detuning $\Delta\lambda_{jk}$ [pm]')
b.set_ylabel(r'$\delta\lambda_{k\leftarrow j}$ [pm]')
panel_title(b, 'b', 'Rule A as a placement criterion')
b.legend(loc='upper right', bbox_to_anchor=(0.99, 0.58), fontsize=5.0,
         handlelength=1.6, labelspacing=0.18, borderaxespad=0.0)

ins = b.inset_axes([0.46, 0.05, 0.52, 0.46])
ins.set_xlim(0, 10)
ins.set_ylim(0, 5)
ins.axis('off')
ins.plot([0.2, 9.8], [2.5, 2.5], color='black', lw=3.0, solid_capstyle='butt',
         zorder=1)
EX = [(1.0, FS.VERM, '$j{=}k{-}2$', '+150 pm', 'forbidden'),
      (4.0, FS.GREEN, '$j{=}k{-}1$', '0', 'allowed'),
      (6.1, FS.ORANGE, '$k$', 'read', ''),
      (8.5, FS.GREEN, '$k{+}1$', 'behind', 'allowed')]
for x0, colr, lab, det, verdict in EX:
    ins.add_patch(Rectangle((x0 - 0.42, 1.75), 0.84, 1.5, facecolor=colr,
                            edgecolor='none', zorder=2))
    ins.text(x0, 2.5, lab, ha='center', va='center', fontsize=4.3,
             color='white', zorder=3, rotation=90)
    ins.text(x0, 1.35, det, ha='center', va='top', fontsize=4.6,
             color=('#B87F00' if colr == FS.ORANGE else colr))
    if verdict:
        ins.text(x0, 3.55, verdict, ha='center', va='bottom', fontsize=4.2,
                 color=colr)
ins.annotate('', xy=(9.6, 0.5), xytext=(0.4, 0.5),
             arrowprops=dict(arrowstyle='-|>', lw=0.7, color='0.4',
                             mutation_scale=6))
ins.text(5.0, 0.1, 'from the laser', ha='center', va='top', fontsize=4.6,
         color='0.4')
ins.text(5.0, 4.9, r'example: detuning from $k$',
         ha='center', va='top', fontsize=4.8, color='0.25')

# ---------------------------------------------------------------------------
# (c) Rule B against the model: layouts, their mean, the closed form, and the
# mean after two references. Data computed by s23_theory.py (s23_ruleB.npz).
# ---------------------------------------------------------------------------
dat = np.load('figs/s23_ruleB.npz')
c.plot(dat['pos'], dat['err'], '.', ms=2, color='0.75', alpha=0.5)
c.plot(dat['cent'], dat['binned'], 'o', color=FS.BLUE, ms=4.5,
       label='mean, no refs')
c.plot(dat['nu_fine'], dat['lawB_fine'], '-', color=FS.VERM, lw=1.4,
       label='Rule B')
c.plot(dat['cent'], dat['binned_r'], 's-', color=FS.GREEN, ms=3.4, lw=1.0,
       mfc='none', label='mean, 2 refs')
c.set_xlabel(r'position in the band $\nu_k$ [pm]')
c.set_ylabel(r'mean bias $\overline{\delta\lambda}(\nu_k)$ [pm]')
panel_title(c, 'c', 'Rule B and two references')
c.legend(fontsize=5.2, loc='upper left', handlelength=1.6, labelspacing=0.18,
         borderaxespad=0.3)

fig.savefig('figs/fig_s26_laws_concept.pdf', bbox_inches='tight',
            pad_inches=0.025)
fig.savefig('figs/fig_s26_laws_concept.png', dpi=300, bbox_inches='tight',
            pad_inches=0.025)
plt.close(fig)

print('Rule A roots at 1 pm: Delta_lo=%.1f pm, Delta_hi=%.1f pm' % (dlo, dhi))
print('Rule A at least half its maximum between %.0f and %.0f pm '
      '(%.2f and %.2f sigma)' % (hlo, hhi, hlo / SIG, hhi / SIG))
print('saved figs/fig_s26_laws_concept.pdf and .png')
