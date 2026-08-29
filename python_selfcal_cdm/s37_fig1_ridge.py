"""s37_fig1_ridge.py - the wide panel of Fig. 1, panels (c) and (d).

Correction of 2026-08-30. Every ridge in panel (c) was sliced off flat at
the top. The cause was the drawing order, not the data. A ridgeline hides
what lies behind it with an opaque fill, so a peak that rises into the row
above survives only if that row was drawn first. The old panel went from
the bottom up, so each row covered the peak below it and every ridge came
out as a trapezoid cut at exactly one row pitch.

Rows now go from the far side towards the reader, high wavelength first,
and the tallest ridge is a fixed 2.6 row pitches. The overlap then reads as
depth rather than as damage. A ridgeline stays the right chart here, since
the sentence the panel has to carry is one ridge, one grating. A 3-D
surface would hide the small ridges behind the large ones and would put the
delay axis in perspective, and that is the axis a reader has to compare
across columns.

Panel (d) keeps its content: the column of samples under the dashed line
in (c), read out as the reflectance of one grating, which is the quantity
that moves with temperature and strain.

Output: fig1_panels/panel_hero.pdf at 361.7 x 135.8 pt, the size the TikZ
composition in Draft_v2/fig1_concept.tex was laid out around.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import common as C
import figstyle as FS

FS.apply(6.4)

# the grating colours of fig1_concept.tex, so the spine and the ridges agree
VERM = '#C74E0A'
ORAN = '#E8A200'
GREE = '#009E73'
PURP = '#CC79A7'
COLS = [VERM, ORAN, GREE, PURP]

NCH = 127                       # chips, so the delay axis is one code period
DELAYS = [18.0, 47.0, 76.0, 110.0]
DETUNE = [-20.0, 130.0, -60.0, 185.0]     # pm from band centre
M = 25                                    # wavelength steps
FOCUS = 1                                 # the grating panel (d) reads out

lam = np.linspace(-0.5, 0.5, M)           # nm
nu = lam * 1000.0 * C.GHZ_PER_PM          # GHz from band centre
PITCH = (lam[-1] - lam[0]) / (M - 1)
PEAK = 2.6 * PITCH                        # tallest ridge, in axis units
HALF = 7.0                                # chips coloured around each delay
TOP = 1.38 * PEAK                         # headroom above the last row

tau = np.linspace(0, NCH, 900)
rng = np.random.default_rng(7)


def chip_peak(t, t0):
    """Triangular correlation peak, two chips wide at the base."""
    return np.clip(1.0 - np.abs(t - t0), 0.0, None)


def wiggle(n):
    """Baseline noise, smoothed to the bandwidth of the correlator.

    One independent sample per plotted point would draw a hairy line at this
    print size, and it would also be wrong: the correlator output is band
    limited by the chip rate, so the noise has to vary on the scale of a
    chip, not of a pixel.
    """
    ker = np.exp(-0.5 * (np.arange(-18, 19) / 6.0) ** 2)
    ker /= ker.sum()
    return np.convolve(rng.normal(size=n + 36), ker, mode='same')[18:18 + n]


def amp(m, k):
    return C.fbg_gauss(nu[m], DETUNE[k] * C.GHZ_PER_PM, C.FBG_FWHM_GHZ)


fig = plt.figure(figsize=(5.023, 1.886))
axc = fig.add_axes([0.072, 0.150, 0.686, 0.800])
axd = fig.add_axes([0.861, 0.150, 0.125, 0.800])

# ------------------------------- panel (c) --------------------------------
# far side first, so a ridge may rise into the rows behind it and stay whole
for m in range(M - 1, -1, -1):
    base = lam[m]
    trace = np.zeros_like(tau)
    for k in range(len(DELAYS)):
        trace += amp(m, k) * chip_peak(tau, DELAYS[k])
    trace = trace + 0.075 * wiggle(tau.size)
    y = base + PEAK * trace
    z = 10 + 2.0 * (M - m)
    axc.fill_between(tau, base - 0.35 * PITCH, y, color='white', lw=0,
                     zorder=z)
    axc.plot(tau, y, color='0.55', lw=0.4, zorder=z + 0.2)
    for k, (d, col) in enumerate(zip(DELAYS, COLS)):
        a = amp(m, k)
        if a < 0.05:
            continue
        w = np.abs(tau - d) <= HALF
        axc.fill_between(tau[w], base, y[w], color=col, alpha=0.32, lw=0,
                         zorder=z + 0.4)
        axc.plot(tau[w], y[w], color=col, lw=0.9, zorder=z + 0.6)

axc.axvline(DELAYS[FOCUS], color='0.15', ls=(0, (3.5, 2.5)), lw=0.8,
            zorder=200)
axc.text(NCH * 0.99, lam[-1] + TOP * 0.97, 'one ridge $=$ one grating',
         fontsize=6.4, color='0.25', ha='right', va='top', zorder=210)

YLO = lam[0] - 1.2 * PITCH
YHI = lam[-1] + TOP
axc.set_xlim(0, NCH)
axc.set_ylim(YLO, YHI)
axc.set_xticks(DELAYS)
axc.set_xticklabels([r'$\tau_%d$' % (i + 1) for i in range(len(DELAYS))])
for lbl, col in zip(axc.get_xticklabels(), COLS):
    lbl.set_color(col)
axc.set_yticks([-0.3, 0.0, 0.3])
axc.set_xlabel(r'delay $\tau$ (chips)', labelpad=1.5)
axc.set_ylabel(r'wavelength $\lambda_m$ (nm)', labelpad=1.5)
axc.tick_params(length=2.0, pad=1.6)

# --------------------------- the cut, in the gutter -----------------------
lam_b = DETUNE[FOCUS] / 1000.0
arr_y = 0.150 + 0.800 * (lam_b - YLO) / (YHI - YLO)
fig.text(0.808, 0.860, '4', fontsize=5.5, color='0.15', ha='center',
         va='center', zorder=6,
         bbox=dict(boxstyle='circle,pad=0.42', facecolor='white',
                   edgecolor='0.35', linewidth=0.5))
fig.text(0.808, 0.750, 'cut', fontsize=6.2, color='0.35', ha='center',
         va='center')
fig.text(0.808, arr_y - 0.105, r'$\tau_2$', fontsize=6.6, color=ORAN,
         ha='center', va='center')
fig.add_artist(matplotlib.patches.FancyArrowPatch(
    (0.774, arr_y), (0.850, arr_y), transform=fig.transFigure,
    arrowstyle='-|>', mutation_scale=5.5, lw=0.8, color='0.35',
    shrinkA=0, shrinkB=0))

# ------------------------------- panel (d) --------------------------------
lam_fine = np.linspace(lam[0], lam[-1], 400)
r_fine = C.fbg_gauss(lam_fine * 1000.0 * C.GHZ_PER_PM,
                     DETUNE[FOCUS] * C.GHZ_PER_PM, C.FBG_FWHM_GHZ)
r_samp = np.array([amp(m, FOCUS) for m in range(M)])

axd.plot(r_fine, lam_fine, color=ORAN, lw=1.1)
axd.plot(r_samp, lam, 'o', color=ORAN, ms=2.0, mew=0)
axd.axhline(lam_b, color='0.15', ls=(0, (3.5, 2.5)), lw=0.8)
axd.text(0.05, lam_b + 0.030, r'$\lambda_B$', fontsize=6.6, color='0.15',
         ha='left', va='bottom')
axd.text(0.5, lam[-1] + 0.5 * TOP,
         r'$\Delta\lambda_B \propto \Delta T,\ \epsilon$',
         fontsize=6.2, color='0.25', ha='center', va='center')

axd.set_xlim(-0.10, 1.16)
axd.set_ylim(YLO, YHI)
axd.set_xticks([0, 1])
axd.set_yticks([])
axd.spines['left'].set_visible(False)
axd.set_xlabel('reflectivity', labelpad=1.5)
axd.tick_params(length=2.0, pad=1.6)

os.makedirs('figs/fig1_panels', exist_ok=True)
fig.savefig('figs/fig1_panels/panel_hero.pdf', pad_inches=0.0,
            transparent=True)
fig.savefig('figs/fig1_panels/panel_hero.png', dpi=300, pad_inches=0.0)
plt.close(fig)

# ridges sit in separate delay bins, so the tallest one is a max, not a sum
tallest = max(amp(m, k) for m in range(M) for k in range(len(DELAYS)))
top_row = max(amp(M - 1, k) for k in range(len(DELAYS)))
print('%d rows, pitch %.4f nm, tallest ridge %.2f pitches'
      % (M, PITCH, PEAK * tallest / PITCH))
print('clearance above the last row: %.2f pitches'
      % ((TOP - PEAK * top_row) / PITCH))
print('saved figs/fig1_panels/panel_hero.pdf and .png')
