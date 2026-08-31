"""s47_overlap_inset.py - the crowded-spectra inset of the measurement-chain
figure, replacing the inset_overlap of s24_concept v11 (commit 354aefb).

Four grating lines R_k(lambda) of equal width crowd one wavelength band, at
the same centres as before (-9, +3, -3, +7 GHz around the band centre), in
the grating colours used by the fiber sketch and by the echo inset. New in
this version: the two axes are drawn and named (wavelength along x, reflected
power along y), every line is labelled FBG1 .. FBG4, and the wavelength the
laser probes in the current step, lambda_m, is marked in the laser colour
with a dot on every line. Those dot heights are the echo amplitudes
A_k(lambda_m) of the echo inset, which is what ties the two panels together.

Output: figs/fig1_panels/inset_overlap.pdf (and .png, .svg), consumed by
Draft_v2/fig1_concept.tex at a width of 5.0 cm. Drawn at that width.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import common as C
import figstyle as FS

FS.apply()

ECOL = ['#C74E0A', '#E8A200', '#009E73', '#CC79A7']
VOLT = '#7d3c98'                            # laser wavelength, as in s35
SIG = C.FBG_FWHM_GHZ / 2.35482
CENT = [-9.0, 3.0, -3.0, 7.0]               # GHz, crowded around one lambda_B
LAM_M = -0.5                                # the probed wavelength, this step
FS_LBL = 5.8
FS_SML = 5.0

lam = np.linspace(-40, 40, 600)
fig = plt.figure(figsize=(1.97, 1.15))      # 5.0 cm wide
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])

# label heights staggered so that four labels over four nearby peaks stay apart
lab_y = {0: 1.08, 2: 1.22, 1: 1.36, 3: 1.50}
for k, (c0, col) in enumerate(zip(CENT, ECOL)):
    g = np.exp(-0.5 * ((lam - c0) / SIG) ** 2)
    ax.plot(lam, g, color=col, lw=0.8)
    ax.fill_between(lam, g, color=col, alpha=0.07, lw=0)
    ax.plot([c0, c0], [1.0, lab_y[k] - 0.03], color=col, lw=0.4, alpha=0.7)
    ax.text(c0, lab_y[k], 'FBG%d' % (k + 1), fontsize=FS_SML, color=col,
            ha='center', va='bottom')
    a = np.exp(-0.5 * ((LAM_M - c0) / SIG) ** 2)
    ax.plot(LAM_M, a, 'o', color=col, ms=2.6, mec='white', mew=0.35, zorder=5)

# the probed wavelength, laser colour
ax.axvline(LAM_M, color=VOLT, lw=0.6, ls=(0, (2.5, 1.8)), ymin=0.0, ymax=0.62)
ax.text(LAM_M + 1.2, -0.05, r'$\lambda_m$', fontsize=FS_LBL, color=VOLT,
        ha='left', va='top')
ax.text(LAM_M + 1.2, -0.30, 'probed in step $m$', fontsize=FS_SML, color=VOLT,
        ha='left', va='top')
ax.text(-38.5, 0.86, 'dots: $A_k(\\lambda_m)$', fontsize=FS_SML, color='0.35',
        ha='left', va='bottom')

# axes: wavelength along x, reflected power along y
X0, Y0 = -40.0, 0.0
ax.annotate('', xy=(41.5, Y0), xytext=(X0, Y0),
            arrowprops=dict(arrowstyle='-|>', lw=0.7, color='0.3',
                            mutation_scale=6))
ax.text(42.0, Y0, r'$\lambda$', fontsize=FS_LBL, color='0.3', ha='left',
        va='center')
ax.annotate('', xy=(X0, 1.62), xytext=(X0, Y0),
            arrowprops=dict(arrowstyle='-|>', lw=0.7, color='0.3',
                            mutation_scale=6))
ax.text(X0 - 0.5, 1.64, r'$P$ (reflected power)', fontsize=FS_LBL,
        color='0.3', ha='left', va='bottom')

ax.set_xlim(X0 - 1.0, 47.0)
ax.set_ylim(-0.62, 1.90)
ax.axis('off')

os.makedirs('figs/fig1_panels', exist_ok=True)
fig.savefig('figs/fig1_panels/inset_overlap.pdf', bbox_inches='tight',
            pad_inches=0.01, transparent=True)
fig.savefig('figs/fig1_panels/inset_overlap.png', dpi=400, bbox_inches='tight',
            pad_inches=0.01)
fig.savefig('figs/fig1_panels/inset_overlap.svg', bbox_inches='tight',
            pad_inches=0.01, transparent=True)
plt.close(fig)
print('A_k(lambda_m):', np.round(np.exp(-0.5 * ((LAM_M - np.array(CENT)) / SIG) ** 2), 3))
print('saved figs/fig1_panels/inset_overlap.pdf, .png, .svg')
