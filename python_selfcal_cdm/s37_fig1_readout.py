"""s37_fig1_readout.py - the lower half of Fig. 1, panels (b), (c) and (d).

Rewritten on 2026-08-30, replacing the ridgeline version. The old panel (c)
drew all 25 correlated rows as a ridge plot. It was large, it repeated the
same shape 25 times, and that repetition carried nothing a reader did not
already have after the first row. What the figure has to say is narrower:
correlate one row, read the height at one delay bin, and that height is one
point of one grating's spectrum. Panel (c) therefore shows a single
correlation, and panel (d) shows the point it contributes with the rest of
the spectrum already filled in.

The three panels come from one simulated record, which is why they are
generated together. The raw stamp in (b) is the record that (c) correlates,
so a reader comparing them is comparing the same numbers.

Outputs, both consumed by Draft_v2/fig1_concept.tex:
  fig1_panels/panel_stamp.pdf  the raw record, one pixel per chip
  fig1_panels/panel_hero.pdf   panels (c) and (d) side by side
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import common as C
import figstyle as FS

FS.apply(6.4)

# the grating colours of fig1_concept.tex, so the spine and the peaks agree
VERM = '#D55E00'
ORAN = '#E69F00'
GREE = '#009E73'
PURP = '#CC79A7'
COLS = [VERM, ORAN, GREE, PURP]

NCH = 127                                 # chips, one code period
DELAYS = [18, 47, 76, 110]                # delay bins, that is grating order
DETUNE = [-20.0, 130.0, -60.0, 185.0]     # pm from band centre
R = 0.10
M = 25                                    # wavelength steps
FOCUS = 1                                 # the grating (d) reads out
SHOWN = 12                                # the step (c) correlates,
#                                           band centre, so all four gratings
#                                           still reflect something

lam = np.linspace(-0.5, 0.5, M)           # nm
nu = lam * 1000.0 * C.GHZ_PER_PM


def amp(m, k):
    return R * C.fbg_gauss(nu[m], DETUNE[k] * C.GHZ_PER_PM, C.FBG_FWHM_GHZ)


# --- the record, and the correlation that reads it ------------------------
code01 = C._mls01(7).astype(float)     # unipolar drive, what the laser does
replica = 2.0 * code01 - 1.0           # bipolar replica, an on chip is +1
rng = np.random.default_rng(4)

record = np.zeros((M, NCH))
for m in range(M):
    for k, d in enumerate(DELAYS):
        record[m] += amp(m, k) * np.roll(code01, d)
record += 0.010 * rng.normal(size=record.shape)

# mean removed before correlating, as in the baseline receiver of Section III-C
corr = np.array([C.periodic_xcorr(record[m] - record[m].mean(), replica)
                 for m in range(M)])
corr *= 2.0                            # unipolar drive halves the peak

# ==========================================================================
# panel (b): the raw record
# ==========================================================================
figb, axb = plt.subplots(figsize=(1.34, 1.10))
axb.imshow(record, aspect='auto', cmap='Greys', interpolation='nearest',
           origin='lower', extent=(0, NCH, 0, M),
           vmin=0.0, vmax=np.percentile(record, 98))
axb.set_xticks([]); axb.set_yticks([])
for sp in axb.spines.values():
    sp.set_visible(True); sp.set_linewidth(0.5); sp.set_color('0.45')
axb.set_xlabel('chip', fontsize=5.6, labelpad=1.2)
axb.set_ylabel(r'step $\lambda_m$', fontsize=5.6, labelpad=1.2)

os.makedirs('figs/fig1_panels', exist_ok=True)
figb.savefig('figs/fig1_panels/panel_stamp.pdf', bbox_inches='tight',
             pad_inches=0.01)
figb.savefig('figs/fig1_panels/panel_stamp.png', dpi=600, bbox_inches='tight',
             pad_inches=0.01)
plt.close(figb)

# ==========================================================================
# panels (c) and (d)
# ==========================================================================
# panele zajmuja dolne dwie trzecie, gora zostaje na luk laczacy
fig = plt.figure(figsize=(4.30, 1.55))
axc = fig.add_axes([0.088, 0.225, 0.505, 0.565])
axd = fig.add_axes([0.735, 0.225, 0.250, 0.565])

# --- (c) one row, correlated ----------------------------------------------
tau = np.arange(NCH)
row = corr[SHOWN]
axc.plot(tau, row, color='0.45', lw=0.6)
for k, (d, col) in enumerate(zip(DELAYS, COLS)):
    w = slice(max(d - 4, 0), d + 5)
    axc.plot(tau[w], row[w], color=col, lw=1.3)
peak = row[DELAYS[FOCUS]]
axc.plot([DELAYS[FOCUS], DELAYS[FOCUS]], [0, peak], color=ORAN, lw=0.6,
         ls=(0, (2.5, 2)), zorder=2)
axc.plot([DELAYS[FOCUS]], [peak], 'o', color=ORAN, ms=3.4, mec='white',
         mew=0.5, zorder=6)

axc.set_xlim(-3, NCH + 2)
axc.set_ylim(-0.012, 1.32 * row.max())
axc.set_xticks(DELAYS)
axc.set_xticklabels([r'$\tau_%d$' % (i + 1) for i in range(len(DELAYS))])
for lbl, col in zip(axc.get_xticklabels(), COLS):
    lbl.set_color(col)
axc.set_yticks([0.0, 0.05, 0.10])
axc.set_xlabel(r'delay $\tau$ [chips]', labelpad=1.2)
axc.set_ylabel('reflectance', labelpad=1.2)
axc.tick_params(length=2.0, pad=1.4)
axc.text(0, 1.28 * row.max(), r'one step, $\lambda_m$', fontsize=6.0,
         color='0.35', ha='left', va='top')

# --- (d) the spectrum that delay bin builds -------------------------------
lam_fine = np.linspace(lam[0], lam[-1], 400)
r_fine = R * C.fbg_gauss(lam_fine * 1000.0 * C.GHZ_PER_PM,
                         DETUNE[FOCUS] * C.GHZ_PER_PM, C.FBG_FWHM_GHZ)
r_samp = corr[:, DELAYS[FOCUS]]

axd.plot(lam_fine, r_fine, color=ORAN, lw=1.0)
axd.plot(lam, r_samp, 'o', color=ORAN, ms=1.9, mew=0)
axd.plot([lam[SHOWN]], [r_samp[SHOWN]], 'o', color=ORAN, ms=3.4, mec='white',
         mew=0.5, zorder=6)
lam_b = DETUNE[FOCUS] / 1000.0
axd.axvline(lam_b, color='0.25', ls=(0, (2.5, 2)), lw=0.7)
axd.text(lam_b + 0.05, 1.17 * R, r'$\lambda_B$', fontsize=6.4, color='0.15',
         ha='left', va='center')
axd.text(lam[0] - 0.02, 0.30 * R, r'$\Delta\lambda_B$' '\n'
         r'$\propto \Delta T,\ \epsilon$', fontsize=6.0,
         color='0.35', ha='left', va='center')

axd.set_xlim(lam[0] - 0.05, lam[-1] + 0.05)
axd.set_ylim(-0.012, 1.32 * R)
axd.set_xticks([-0.4, 0.0, 0.4])
axd.set_yticks([])
axd.spines['left'].set_visible(False)
axd.set_xlabel(r'$\lambda_m$ [nm]', labelpad=1.2)
axd.tick_params(length=2.0, pad=1.4)

# the one link the figure is about: this peak is that point
fig.add_artist(ConnectionPatch(
    xyA=(DELAYS[FOCUS], peak), coordsA=axc.transData,
    xyB=(lam[SHOWN], r_samp[SHOWN]), coordsB=axd.transData,
    arrowstyle='-|>', mutation_scale=6, lw=0.8, color='0.35',
    connectionstyle='arc3,rad=-0.38', shrinkA=3.0, shrinkB=2.0))
fig.text(0.664, 0.40, 'repeat over' '\n' r'every $\lambda_m$', fontsize=5.8,
         color='0.35', ha='center', va='center')

fig.savefig('figs/fig1_panels/panel_hero.pdf', bbox_inches='tight',
            pad_inches=0.005, transparent=True)
fig.savefig('figs/fig1_panels/panel_hero.png', dpi=300, bbox_inches='tight',
            pad_inches=0.005)
plt.close(fig)

print('peak at bin %d reads %.4f, the model puts %.4f there'
      % (DELAYS[FOCUS], peak, amp(SHOWN, FOCUS)))
print('spectrum top %.4f against R = %.2f' % (r_samp.max(), R))
print('saved figs/fig1_panels/panel_stamp.pdf and panel_hero.pdf')
