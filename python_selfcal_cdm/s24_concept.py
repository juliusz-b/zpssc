"""s24_concept.py - the whole system in one figure, JLT-Fig.1 style.

  a  the hardware chain with micro-insets of the real signals at the points
     where they live: the code at the source, the delayed echoes on the fiber,
     the sum at the photodiode. Small, computed, attached by thin leaders.
  b  why wavelength alone cannot separate the gratings: their spectra share one
     band and overlap; the delay axis is what tells them apart.
  c  what the correlator returns: the delay-wavelength map (full model).
  d  the readout: one column is one grating's spectrum, its peak is lambda_B.

Colour continuity carries the argument: each grating keeps its colour from the
fiber, through the overlap panel and the map, to the recovered spectrum.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS

FS.apply()

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
NCH = 127
_MS = 1.0 - 2.0 * C._mls01(7)
_MS = np.roll(_MS, -int(np.argmax(_MS == 1.0)))
ACORR = C.periodic_xcorr(_MS, _MS)
ECOL = [FS.VERM, FS.ORANGE, FS.GREEN, FS.PURPLE]

# --- signals for the micro-insets -------------------------------------------
OSR = 8
D01 = 0.5 * (1 + _MS)
drv = np.tile(np.repeat(D01, OSR), 2)
tt = np.arange(len(drv)) / OSR
taus = [4.0, 9.0, 14.0]
amps = [1.0, 0.72, 0.5]
ech = []
for tc, a in zip(taus, amps):
    e = np.zeros_like(drv); sh = int(tc * OSR)
    e[sh:] = a * drv[:len(drv) - sh]
    ech.append(e)
rngN = np.random.default_rng(3)
sm = np.sum(ech, axis=0) + rngN.normal(0, 0.05, len(drv))
WIN = tt < 22

# --- map data (full model: shadowing + leakage) ------------------------------
M = 120
nu = np.linspace(-2.3 * F, 2.3 * F, M)
K = 4
bins = np.array([16, 47, 76, 105])
nub = np.array([-16.0, 9.0, -5.0, 19.0])
R = 0.05
shapes = np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
tcum = np.ones((K, M))
for k in range(1, K):
    tcum[k] = tcum[k - 1] * (1.0 - R * shapes[k - 1]) ** 2
prim = R * shapes * tcum
grid = np.zeros((NCH, M))
grid[bins] = prim
Wl = ACORR[(np.arange(NCH)[:, None] - bins[None, :]) % NCH]
grid = grid + Wl @ prim
HL = 1
spec = grid[bins[HL]]
mu_fit = C.gauss_fit_peak(nu, spec)

# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(7.1, 2.90))
bg = fig.add_axes([0, 0, 1, 1]); bg.set_xlim(0, 100); bg.set_ylim(0, 100)
FS.despine_all(bg)
TXT = '0.25'
GY = '0.45'


def mini(rect):
    a = fig.add_axes(rect)
    for sp in a.spines.values():
        sp.set_visible(False)
    a.set_xticks([]); a.set_yticks([])
    return a


# =========================== a: the chain ===================================
YW = 76
# source box with the sweep glyph inside
bg.add_patch(FancyBboxPatch((2.0, YW - 10), 9.4, 20, fc='white', ec='0.35',
                            lw=0.7, boxstyle='round,pad=0.3,rounding_size=0.9'))
bg.plot([3.4, 5.9, 5.9, 8.4], [YW + 3.2, YW + 7.4, YW + 3.2, YW + 7.4],
        color=FS.BLUE, lw=0.8)
bg.text(9.6, YW + 5.3, r'$\lambda(t)$', fontsize=5.2, color=FS.BLUE,
        va='center')
bg.text(6.7, YW - 4.5, 'swept\nVCSEL', fontsize=5.8, ha='center', va='center',
        color='0.15')
# circulator
bg.plot([11.6, 14.3], [YW, YW], color='0.35', lw=0.8)
bg.plot([15.7], [YW], marker='o', ms=8.0, mfc='white', mec='0.35', mew=0.7)
# fiber, gratings (labels BELOW the fiber), dashed continuation
bg.plot([17.1, 55.0], [YW, YW], color='0.35', lw=1.1, solid_capstyle='round')
bg.plot([55.0, 60.5], [YW, YW], color='0.35', lw=1.1, ls=(0, (1.2, 1.6)))
gx = [24.0, 31.0, 38.0, 52.0]
for x, col, i in zip(gx, ECOL, range(4)):
    for dx in (-0.35, 0.0, 0.35):
        bg.plot([x + dx, x + dx], [YW - 3.2, YW + 3.2], color=col, lw=1.0)
    bg.text(x, YW - 7.8, '$z_%d$' % (i + 1), fontsize=5.6, ha='center',
            va='top', color=col)
bg.text(45.0, YW - 7.8, r'$\cdots$', fontsize=7, ha='center', va='top',
        color='0.35')
bg.text(62.0, YW, 'FBG array\nin one band', fontsize=5.8, color=TXT,
        va='center')
# photodiode (label on the left) and correlator
bg.plot([15.7, 15.7], [YW - 3.2, YW - 12], color='0.35', lw=0.8)
bg.plot([15.7], [YW - 14.5], marker='v', ms=4.2, color='0.35')
bg.plot([14.6, 16.8], [YW - 17.6, YW - 17.6], color='0.35', lw=0.8)
bg.text(13.4, YW - 15.0, 'PD', fontsize=5.6, color=TXT, ha='right',
        va='center')
bg.plot([15.7, 15.7], [YW - 19.4, YW - 24.5], color='0.35', lw=0.8)
bg.plot([15.7], [YW - 28.5], marker='o', ms=10.0, mfc='white', mec='0.35',
        mew=0.7)
bg.text(15.7, YW - 28.9, r'$\otimes c$', fontsize=6.0, ha='center', va='center',
        color='0.15')
bg.text(19.6, YW - 28.8, 'correlate', fontsize=5.4, color=TXT, ha='left',
        va='center')

# micro-inset: the code, above the source wire
axcode = mini([0.135, 0.845, 0.095, 0.115])
axcode.plot(tt[WIN], drv[WIN], color=FS.BLUE, lw=0.65)
bg.text(18.2, 97.0, 'code $c(t)$', fontsize=5.4, color=FS.BLUE, ha='center')
bg.plot([13.6, 13.0], [84.2, YW + 1.6], color='0.65', lw=0.45)

# micro-inset: the three echoes, above the fiber
axech = mini([0.30, 0.845, 0.135, 0.115])
for e, col, off in zip(ech, ECOL, (0.0, 1.25, 2.5)):
    axech.plot(tt[WIN], off + e[WIN], color=col, lw=0.65)
axech.set_ylim(-0.2, 3.8)
bg.text(36.7, 97.0, 'delayed echoes', fontsize=5.4, color=TXT, ha='center')
bg.plot([33.5, 31.5], [84.2, YW + 3.4], color='0.65', lw=0.45)

# micro-inset: the sum at the photodiode, in the free column under the source
axsum = mini([0.015, 0.50, 0.105, 0.10])
axsum.plot(tt[WIN], sm[WIN], color='0.3', lw=0.6)
bg.text(6.6, 61.8, 'sum at PD', fontsize=5.4, color=TXT, ha='center')
bg.plot([12.3, 14.4], [56.0, YW - 14.5], color='0.65', lw=0.45)

# arrow: correlator -> map panel
bg.annotate('', xy=(30.0, 42.5), xytext=(17.6, YW - 32.5),
            arrowprops=dict(arrowstyle='-|>', lw=1.0, color=GY,
                            connectionstyle='arc3,rad=-0.18'))

# =========================== b: overlapping spectra ==========================
axO = fig.add_axes([0.05, 0.145, 0.15, 0.26])
lamg = np.linspace(-70, 70, 300)
for c0, col in zip(nub, ECOL):
    g = np.exp(-0.5 * ((lamg - c0) / SIG) ** 2)
    axO.plot(lamg * PM / 1000, g, color=col, lw=0.8)
    axO.fill_between(lamg * PM / 1000, g, color=col, alpha=0.12, lw=0)
axO.set_xlim(-0.56, 0.56); axO.set_ylim(0, 1.62)
axO.set_xticks([-0.4, 0, 0.4]); axO.set_yticks([])
axO.tick_params(labelsize=5.4, pad=1.5)
axO.set_xlabel(r'$\lambda$ offset (nm)', fontsize=5.8, labelpad=1)
for sp in ('top', 'right', 'left'):
    axO.spines[sp].set_visible(False)
axO.text(0.5, 0.97, 'spectra overlap: wavelength\nalone cannot separate them',
         transform=axO.transAxes, fontsize=5.6, ha='center', va='top',
         color='0.15')

# =========================== c: the map =====================================
axB = fig.add_axes([0.315, 0.145, 0.355, 0.335])
axB.imshow(grid.T, aspect='auto', origin='lower', cmap='Greys',
           extent=[0, NCH, nu[0] * PM / 1000, nu[-1] * PM / 1000],
           vmax=grid.max() * 0.9)
for b, col in zip(bins, ECOL):
    axB.plot(b, 0.505, marker='v', ms=3.0, color=col, clip_on=False)
axB.axvspan(bins[HL] - 3.5, bins[HL] + 3.5, color=FS.ORANGE, alpha=0.16, lw=0)
axB.set_xlabel(r'delay $\tau$  $\rightarrow$  which grating', fontsize=6.0,
               labelpad=1.5)
axB.set_ylabel(r'$\lambda$ offset (nm)', fontsize=6.0, labelpad=0.5)
axB.tick_params(labelsize=5.4, pad=1.5)
axB.set_yticks([-0.4, 0, 0.4]); axB.set_xticks([0, 40, 80, 120])
for sp in ('top', 'right'):
    axB.spines[sp].set_visible(False)

# arrow map -> readout
bg.annotate('', xy=(73.0, 30), xytext=(68.4, 30),
            arrowprops=dict(arrowstyle='-|>', lw=1.0, color=GY))
bg.text(70.7, 36.5, 'one\ncolumn', fontsize=5.3, ha='center', color=TXT)

# =========================== d: the readout =================================
axC = fig.add_axes([0.80, 0.145, 0.185, 0.335])
lam = nu * PM / 1000
axC.plot(lam, spec / spec.max(), color=ECOL[HL], lw=1.0)
axC.axvline(mu_fit * PM / 1000, color='0.35', ls=(0, (3, 2)), lw=0.7)
axC.text(mu_fit * PM / 1000 + 0.06, 0.99, r'$\lambda_B$', fontsize=6.6,
         color='0.15', va='top')
axC.text(0.03, 0.86, r'$\Delta\lambda_B \propto \Delta T,\ \varepsilon$',
         fontsize=5.8, color='0.15', transform=axC.transAxes)
axC.set_xlabel(r'$\lambda$ offset (nm)', fontsize=6.0, labelpad=1.5)
axC.set_ylabel('refl.', fontsize=6.0, labelpad=0.5)
axC.set_xlim(-0.55, 0.55); axC.set_ylim(0, 1.06)
axC.set_xticks([-0.4, 0, 0.4]); axC.set_yticks([0, 1])
axC.tick_params(labelsize=5.4, pad=1.5)

FS.panel(fig, 0.006, 0.985, 'a')
FS.panel(fig, 0.006, 0.47, 'b')
FS.panel(fig, 0.272, 0.53, 'c')
FS.panel(fig, 0.752, 0.53, 'd')

fig.savefig('figs/fig_s24_concept.png', bbox_inches='tight', pad_inches=0.02)
print('fitted lambda_B of the highlighted grating: %.1f pm (true %.1f pm)'
      % (mu_fit * PM, nub[HL] * PM))
print('saved figs/fig_s24_concept.png')
