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
    FS.letter(ax, letter)


def gaussian(x, amp, centre, width, baseline):
    return baseline + amp * np.exp(-0.5 * ((x - centre) / width) ** 2)


fig = plt.figure(figsize=(3.5, 1.9))
gs0 = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.40,
                       left=0.10, right=0.995, bottom=0.19, top=0.90)
pa = fig.add_subplot(gs0[0, 0])
b = fig.add_subplot(gs0[0, 1])

# ---------------------------------------------------------------------------
# (a) the three regimes merged on one axis, the light-path sketch as an
# inset in the empty left part (x extended to -730 pm to make that room)
# ---------------------------------------------------------------------------
panel_title(pa, 'a', 'Rule A on the spectrum')
nu = np.linspace(-330, 560, 1800)
wanted = np.exp(-0.5 * (nu / SIG) ** 2)
CASES = ((0.0, FS.BLUE, 'co-tuned', '-'), (DSTAR, FS.VERM, 'worst', (0, (4, 1.5, 1, 1.5))),
         (400.0, FS.GREEN, 'far', (0, (1.2, 1.6))))
ANN = ((-195, 0.70), (300, 0.62), (435, 0.30))
shifts = []
for (det, colr, name, lsty), (tx, ty) in zip(CASES, ANN):
    two_pass = (1.0 - R * np.exp(-0.5 * ((nu - det) / SIG) ** 2)) ** 2
    received = wanted * two_pass
    (amp, mu, sig_, base), _ = curve_fit(
        gaussian, nu, received, p0=[0.9, 0.0, SIG, 0.0], maxfev=20000)
    shifts.append(mu)
    pa.plot(nu, two_pass, color=colr, lw=0.85, alpha=0.85, ls=lsty, zorder=2)
    pa.fill_between(nu, received, wanted, where=wanted >= received,
                    color=colr, alpha=0.16, lw=0, zorder=1)
    pa.plot(nu, received, color=colr, lw=1.4, ls=lsty, zorder=4)
    pa.axvline(mu, color=colr, lw=0.7, ls=(0, (1.5, 2)), zorder=1)
    lab = ('%.1f' % mu).replace('-0.0', '0.0')
    # leader lines to both curves of the case: the thin two-pass transmission and the thick return
    ARROW = {'co-tuned': ((-120.0, 0.90), (-60.0, 0.71)), 'worst': ((130.0, 0.81), (80.0, 0.62)), 'far': ((400.0, 0.81), (120.0, 0.52))}
    for xy in ARROW[name]:
        pa.annotate('', xy=xy, xytext=(tx, ty), arrowprops=dict(arrowstyle='-', lw=0.6, ls=lsty, color=colr, shrinkA=9, shrinkB=1))
    pa.text(tx, ty, name + '\n' + r'$%s$ pm' % lab, fontsize=6, color=colr, ha='center', va='center',
            bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
pa.plot(nu, wanted, color='0.12', ls=(0, (4, 2.5)), lw=1.2, zorder=6)
pa.set_xlim(-330, 560)
pa.set_ylim(0.0, 1.36)
pa.set_yticks([0, 0.5, 1.0])
pa.set_xticks([-200, 0, 200, 400])
pa.set_xlabel(r'wavelength offset from $\lambda_{B,k}$ [pm]')
pa.set_ylabel('normalized reflectance')
roles = [Line2D([], [], color='0.3', lw=0.85, alpha=0.85),
         Line2D([], [], color='0.12', ls=(0, (4, 2.5)), lw=1.2),
         Line2D([], [], color='0.3', lw=1.4)]
pa.legend(roles, [r'$(1-R_j)^2$', r'$R_k(\lambda)$', r'$A_k$'],
          fontsize=5.5, loc='upper center', ncol=3, handlelength=1.8, borderaxespad=0.25,
          handletextpad=0.4, columnspacing=0.8)


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
b.plot(D, bias, color=FS.C_THEORY, lw=1.4, label='shadowing shift')
# the full two-pass model with an unrestricted Gaussian fit, as markers
Dm = np.linspace(25.0, 625.0, 13)
mod = []
for d in Dm:
    rec = wanted * (1.0 - R * np.exp(-0.5 * ((nu - d) / SIG) ** 2)) ** 2
    (amp, mu, sg_, base), _ = curve_fit(
        gaussian, nu, rec, p0=[0.9, 0.0, SIG, 0.0], maxfev=20000)
    mod.append(mu)
b.plot(Dm, mod, 's', ms=3.6, mfc='none', mec=FS.C_MEAS, mew=0.9, ls='none',
       label='full model')
b.axhline(-EPS, color='0.35', lw=0.75, ls=(0, (4, 2)),
          label='tol. $\\epsilon=1$ pm')
b.axvline(dlo, color=FS.VERM, lw=0.8, ls=(0, (2, 2)))
b.axvline(dhi, color=FS.VERM, lw=0.8, ls=(0, (2, 2)))
b.text(dlo + 7, -6.4, r'$\Delta\lambda_{\mathrm{lo}}$', rotation=90, ha='left', va='center', fontsize=6, color=FS.VERM)
b.text(dhi - 7, -2.7, r'$\Delta\lambda_{\mathrm{hi}}$', rotation=90, ha='right', va='center', fontsize=6, color=FS.VERM)
b.plot(DSTAR, bias.min(), 'o', color=FS.VERM, ms=3.5)
b.text(12, -9.75, 'min. $-0.81R\\sigma$',
       fontsize=4.8, color=FS.VERM, va='bottom', ha='left')

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
b.set_ylim(-9.8, 1.15)
b.set_yticks([-8, -4, 0])
b.set_xlabel(r'pair detuning $\Delta\lambda_{jk}$ [pm]')
b.set_ylabel(r'$\delta\lambda_{k\leftarrow j}$ [pm]')
panel_title(b, 'b', 'Rule A as a placement criterion')
b.legend(loc='upper right', bbox_to_anchor=(0.99, 0.585), fontsize=5.0,
         handlelength=1.4, labelspacing=0.15, borderaxespad=0.0, handletextpad=0.4)

ins = b.inset_axes([0.40, 0.03, 0.59, 0.50])
ins.set_xlim(0, 10)
ins.set_ylim(0, 5)
ins.axis('off')
ins.plot([0.2, 9.8], [2.5, 2.5], color='black', lw=3.0, solid_capstyle='butt',
         zorder=1)
EX = [(0.9, FS.VERM, '$j{=}k{-}2$', '+150 pm', 'forbidden'),
      (3.7, FS.GREEN, '$j{=}k{-}1$', '0', 'allowed'),
      (6.4, FS.ORANGE, '$k$', 'read', ''),
      (9.1, FS.GREEN, '$k{+}1$', 'behind', 'allowed')]
for x0, colr, lab, det, verdict in EX:
    ins.add_patch(Rectangle((x0 - 0.42, 1.55), 0.84, 1.9, facecolor=colr,
                            edgecolor='none', zorder=2))
    ins.text(x0, 2.5, lab, ha='center', va='center', fontsize=4.2,
             color='white', zorder=3, rotation=90)
    ins.text(x0, 1.2, det, ha='center', va='top', fontsize=4.2,
             color=('#B87F00' if colr == FS.ORANGE else colr))
ins.annotate('', xy=(9.6, 0.75), xytext=(0.4, 0.75),
             arrowprops=dict(arrowstyle='-|>', lw=0.7, color='0.4',
                             mutation_scale=6))
ins.text(5.0, 0.55, 'from the laser', ha='center', va='top', fontsize=4.4,
         color='0.4')

# Rule B against the model (former panel c) is drawn by s16_principle.py as Fig. 6(b) from figs/s23_ruleB.npz.

fig.savefig('figs/fig_s26_laws_concept.pdf', bbox_inches='tight',
            pad_inches=0.025)
fig.savefig('figs/fig_s26_laws_concept.png', dpi=300, bbox_inches='tight',
            pad_inches=0.025)
plt.close(fig)

print('Rule A roots at 1 pm: Delta_lo=%.1f pm, Delta_hi=%.1f pm' % (dlo, dhi))
print('Rule A at least half its maximum between %.0f and %.0f pm '
      '(%.2f and %.2f sigma)' % (hlo, hhi, hlo / SIG, hhi / SIG))
print('saved figs/fig_s26_laws_concept.pdf and .png')
