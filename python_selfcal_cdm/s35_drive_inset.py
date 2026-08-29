"""s35_drive_inset.py - the drive inset of Fig. 1, with the two control
signals kept apart.

Correction of 2026-08-29. The previous inset drew one staircase with the
code riding on it and labelled the axis I, which says the laser is tuned by
its current. That is not this laser. In a HCG-VCSEL the wavelength is set by
the voltage on the high-contrast-grating mirror, while the drive current is
held at a fixed bias and carries only the intensity modulation. Markowski
et al. state it directly for the same device: the wavelength is set through
the voltage applied to the HCG, fed to a separate input, and the code enters
through an RF port over a bias tee on a constant 12 mA sink.

So the inset now shows two traces. The HCG voltage steps once per wavelength
step and carries no code. The current sits at a constant bias and carries the
code. That separation is also what the rest of the paper assumes: the
wavelength axis is calibrated against a voltage-to-wavelength curve, and the
chirp of Section III-D comes from the current modulation, not from tuning.

Output: fig1_panels/inset_drive.pdf, consumed by Draft_v2/fig1_concept.tex.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import common as C
import figstyle as FS

FS.apply(6.4)

VOLT = '#7d3c98'      # HCG voltage, wavelength control
CURR = '#1f77b4'      # drive current, the code

NSTEP = 4                       # wavelength steps shown
NCHIP = 12                      # chips drawn per step, for legibility
code = C._mls01(4)[:NCHIP]      # a short slice, only to look like a code

t = np.arange(NSTEP * NCHIP)
volt = np.repeat(np.arange(NSTEP), NCHIP).astype(float)
curr = np.tile(code, NSTEP).astype(float)

fig, ax = plt.subplots(2, 1, figsize=(2.05, 0.93), sharex=True,
                       gridspec_kw=dict(height_ratios=[1.0, 1.0], hspace=0.42))

# --- HCG voltage: the wavelength staircase, no code on it -----------------
ax[0].step(t, volt, where='post', color=VOLT, lw=1.5)
ax[0].set_ylim(-0.45, NSTEP - 0.35)
ax[0].text(0.0, NSTEP - 0.42, r'$V_\mathrm{HCG}$: sets $\lambda$',
           fontsize=5.8, color=VOLT, va='top', ha='left')

# --- drive current: constant bias, code on top ----------------------------
ax[1].axhline(0.0, color='0.72', lw=0.6, ls=(0, (2.5, 2)))
ax[1].step(t, curr, where='post', color=CURR, lw=1.1)
ax[1].set_ylim(-0.55, 2.05)
ax[1].text(0.0, 2.0, r'$I$: bias $+$ code', fontsize=5.8, color=CURR,
           va='top', ha='left')
ax[1].text(0.0, -0.52, 'bias', fontsize=5.2, color='0.55', va='bottom')

for a in ax:
    a.set_xlim(0, NSTEP * NCHIP)
    a.set_xticks([]); a.set_yticks([])
    for sp in ('top', 'right', 'left', 'bottom'):
        a.spines[sp].set_visible(False)

os.makedirs('figs/fig1_panels', exist_ok=True)
for target in ('figs/fig1_panels/inset_drive.pdf',):
    fig.savefig(target, bbox_inches='tight', pad_inches=0.01, transparent=True)
fig.savefig('figs/fig1_panels/inset_drive.png', dpi=300, bbox_inches='tight',
            pad_inches=0.01)
plt.close(fig)

print('drive inset: HCG voltage steps the wavelength, current carries the code')
print('saved figs/fig1_panels/inset_drive.pdf and .png')
