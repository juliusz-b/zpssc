"""s24_storyboard.py - the life of one code period, drawn to be read at a glance.

Design brief: the fiber carries a distance axis on top; every signal lane below
shares the matching time axis, t = 2nz/c, so each grating's echo begins directly
beneath the grating and the correlation peak lands on the same vertical. The
three light dashed verticals ARE the method. No decoration, no boxes, no step
bubbles; direct labels, one accent palette, all waveforms computed with the real
m-sequence.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Polygon
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS

FS.apply()

C_LIGHT, N_GROUP = 2.99792458e8, 1.468
B = 25e6
NBITS, N = 7, 127
OSR = 16
MSEQ = 1.0 - 2.0 * C._mls01(NBITS)
MSEQ = np.roll(MSEQ, -int(np.argmax(MSEQ == 1.0)))
DRIVE01 = 0.5 * (1.0 + MSEQ)

z = np.array([16.3, 44.9, 73.4])
refl = np.array([0.10, 0.07, 0.05])
tau = 2.0 * N_GROUP * z / C_LIGHT * B                # 4, 11, 18 chips
CH = 26.0
ECOL = [FS.VERM, FS.ORANGE, FS.GREEN]

drive3 = np.tile(np.repeat(DRIVE01, OSR), 3)
t3 = np.arange(len(drive3)) / OSR
rng = np.random.default_rng(5)
echo3 = []
for r, tc in zip(refl, tau):
    sh = int(round(tc * OSR))
    e = np.zeros_like(drive3)
    e[sh:] = r * drive3[:len(drive3) - sh]
    echo3.append(e)
rx3 = np.sum(echo3, axis=0) + rng.normal(0, 0.005, len(drive3))
rx_ss = rx3[N * OSR:2 * N * OSR].reshape(-1, OSR).mean(axis=1)
corr = C.periodic_xcorr(rx_ss - rx_ss.mean(), MSEQ)
lag = np.arange(N)
sel = t3 < CH

# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.1, 3.55))
ax.set_xlim(-5.6, 28.6)
ax.set_ylim(-0.75, 10.75)
FS.despine_all(ax)

GY = '#B8B8B8'
TXT = '0.25'

# alignment verticals, behind everything
for tc in tau:
    ax.plot([tc, tc], [0.55, 9.30], color=GY, ls=(0, (2, 2)), lw=0.55, zorder=0)

# --- fiber with its distance axis (top) --------------------------------------
YF = 9.65
ax.plot([-1.15, CH + 0.4], [YF, YF], color='0.35', lw=1.3, zorder=2,
        solid_capstyle='round')
# distance ticks on the fiber
for tc, zz in zip(tau, z):
    ax.plot([tc, tc], [YF + 0.10, YF + 0.30], color='0.35', lw=0.7)
    ax.text(tc, YF + 0.90, '$z_%d$ = %.0f m' % (list(tau).index(tc) + 1, zz),
            ha='center', fontsize=6.2, color=TXT)
# FBG symbols: small rectangles with internal grating lines
for tc, col in zip(tau, ECOL):
    ax.add_patch(plt.Rectangle((tc - 0.42, YF - 0.26), 0.84, 0.52, fc='white',
                               ec=col, lw=0.9, zorder=4))
    for dx in (-0.18, 0.0, 0.18):
        ax.plot([tc + dx, tc + dx], [YF - 0.18, YF + 0.18], color=col, lw=0.7,
                zorder=5)
ax.text(CH + 0.55, YF, 'fiber', fontsize=6.2, color=TXT, va='center')

# source and circulator
ax.add_patch(FancyBboxPatch((-5.45, YF - 0.44), 2.85, 0.88,
                            boxstyle='round,pad=0.03,rounding_size=0.10',
                            fc='white', ec='0.35', lw=0.7, zorder=3))
ax.text(-4.02, YF, 'swept VCSEL,
code-modulated', fontsize=5.8, ha='center',
        va='center', color='0.15')
from matplotlib.patches import Arc
ax.add_patch(Circle((-1.85, YF), 0.34, fc='white', ec='0.35', lw=0.7, zorder=4))
ax.add_patch(Arc((-1.85, YF), 0.34, 0.34, theta1=-50, theta2=210, lw=0.55,
                 color='0.35', zorder=5))
ax.add_patch(Polygon([[-1.72, YF - 0.10], [-1.62, YF - 0.20],
                      [-1.77, YF - 0.22]], closed=True, fc='0.35', ec='none',
                     zorder=6))
ax.plot([-2.60, -2.19], [YF, YF], color='0.35', lw=0.8)
# short drop from the circulator to a photodiode, kept clear of the labels
ax.plot([-1.85, -1.85], [YF - 0.34, YF - 0.78], color='0.35', lw=0.8)
ax.add_patch(Polygon([[-2.00, YF - 0.78], [-1.70, YF - 0.78],
                      [-1.85, YF - 1.04]], closed=True, fc='0.35', ec='none',
                     zorder=4))
ax.plot([-2.00, -1.70], [YF - 1.08, YF - 1.08], color='0.35', lw=0.8)
ax.text(-1.38, YF - 0.95, 'PD', fontsize=5.8, color=TXT, va='center')

# --- lane 1: transmitted code ------------------------------------------------
Y1 = 7.95
ax.plot(t3[sel], Y1 + 0.60 * drive3[sel], color=FS.BLUE, lw=0.8)
ax.text(-0.45, Y1 + 0.30, 'transmitted\ncode  $c(t)$', ha='right', fontsize=6.4,
        color='0.15')
i0 = int(np.argmax(drive3 > 0.5))
x0 = t3[i0]
ax.annotate('', xy=(x0 + 1.0, Y1 + 0.86), xytext=(x0, Y1 + 0.86),
            arrowprops=dict(arrowstyle='|-|,widthA=0.12,widthB=0.12', lw=0.55,
                            color=TXT))
ax.text(x0 + 0.5, Y1 + 1.02, r'$T_{\mathrm{chip}}$', fontsize=5.8, ha='center',
        color=TXT)

# --- lanes 2: echoes ---------------------------------------------------------
YE = [6.75, 6.05, 5.35]
for y0, e, tc, col, i in zip(YE, echo3, tau, ECOL, range(3)):
    ax.plot(t3[sel], y0 + 4.6 * e[sel], color=col, lw=0.8)
    ax.plot([tc], [y0 - 0.09], marker='^', ms=2.6, color=col, clip_on=False)
    ax.text(CH + 0.55, y0 + 0.16, r'$\times R_%d$' % (i + 1), fontsize=6.0,
            color=col, va='center')
ax.text(-0.45, YE[1] + 0.35, 'echoes, delayed\nby  $\\tau_k = 2nz_k/c$',
        ha='right', fontsize=6.4, color='0.15')

# --- lane 3: detected sum ----------------------------------------------------
Y3 = 4.15
ax.plot(t3[sel], Y3 + 2.9 * rx3[sel], color='0.30', lw=0.8)
ax.text(-0.45, Y3 + 0.35, 'detected sum\n+ noise', ha='right', fontsize=6.4,
        color='0.15')

# --- matched-filter divider --------------------------------------------------
YD = 3.30
ax.plot([0.0, 9.0], [YD, YD], color=GY, lw=0.6)
ax.plot([17.0, CH], [YD, YD], color=GY, lw=0.6)
ax.text(13.0, YD, r'correlate with  $c(t)$', fontsize=6.2, ha='center',
        va='center', color=TXT, style='italic')

# --- lane 4: correlation -----------------------------------------------------
Y4 = 0.95
m = lag <= CH
ax.plot(lag[m], Y4 + np.maximum(corr[m], -0.002) / corr.max() * 1.55,
        color=FS.BLUE, lw=0.9)
for tc, col in zip(tau, ECOL):
    k = int(round(tc))
    ax.plot(k, Y4 + corr[k] / corr.max() * 1.55 + 0.14, marker='v', ms=2.8,
            color=col)
ax.text(-0.45, Y4 + 0.55, 'correlation\noutput', ha='right', fontsize=6.4,
        color='0.15')
ax.text(CH + 0.55, Y4 + 1.30, 'height =\n$R_k(\\lambda)$ at the\ncurrent sweep '
        '$\\lambda$', fontsize=5.9, color=TXT, va='top')

# bottom delay axis
YA = 0.42
ax.plot([0, CH], [YA, YA], color='0.35', lw=0.7)
for x in range(0, 26, 5):
    ax.plot([x, x], [YA, YA - 0.13], color='0.35', lw=0.6)
    ax.text(x, YA - 0.30, str(x), fontsize=5.8, ha='center', va='top',
            color=TXT)
ax.text(CH / 2, -0.62, r'lag $\tau$ (chips)  $\equiv$  round-trip delay'
        r'  $\equiv$  position', fontsize=6.4, ha='center', color='0.15')

fig.subplots_adjust(left=0.005, right=0.995, top=0.99, bottom=0.01)
fig.savefig('figs/fig_s24_storyboard.png', bbox_inches='tight', pad_inches=0.02)
print('peaks at %s chips, heights %s'
      % ([int(round(x)) for x in tau],
         np.round([corr[int(round(x))] for x in tau], 4).tolist()))
print('saved figs/fig_s24_storyboard.png')
