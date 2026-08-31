"""s35_drive_inset.py - the drive inset of the measurement-chain figure.

Version of 2026-08-31. Three traces in time and one spectrum.

  V_HCG(t)   the voltage on the high-contrast-grating mirror steps once per
             wavelength step. The steps are numbered lambda_1 .. lambda_4, the
             same index m as the rows of the raw record.
  I(t)       the drive current sits at a constant bias and carries the code
             c(t). One code period per wavelength step, the same chips at
             every step. The chips are the first ten of the m-sequence used
             by the echo inset (s24 v11), so the code seen here is the code
             that comes back delayed in the echo inset.
  P_out(t)   the optical output: intensity follows c(t), the wavelength is
             lambda_m during step m.
  spectrum   the optical spectrum of the swept laser: one narrow line per
             step, at lambda_1 .. lambda_4.

In a HCG-VCSEL the wavelength is set by the HCG voltage and the code enters
through the current (Markowski et al.), which is why the two controls are
kept on separate traces.

Output: figs/fig1_panels/inset_drive.pdf (and .png, .svg), consumed by
Draft_v2/fig1_concept.tex at a width of 4.5 cm. The figure is drawn at that
width, so the font sizes below are the printed ones.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import common as C
import figstyle as FS

FS.apply()

VOLT = '#7d3c98'      # HCG voltage, wavelength control
CURR = FS.BLUE        # drive current, the code
LIGHT = '#404040'     # optical output
FS_LBL = 5.8
FS_SML = 5.0

# --- the code: identical construction to the echo inset -------------------
_MS = 1.0 - 2.0 * C._mls01(7)
_MS = np.roll(_MS, -int(np.argmax(_MS == 1.0)))
code01 = 0.5 * (1 + _MS)

NSTEP = 4
NCHIP = 10                              # chips drawn per step
chips = code01[:NCHIP]
print('chips drawn per step:', chips.astype(int))

t = np.arange(NSTEP * NCHIP + 1)
volt = np.repeat(np.arange(NSTEP), NCHIP).astype(float)
volt = np.append(volt, volt[-1])
curr = np.tile(chips, NSTEP)
curr = np.append(curr, curr[-1])
T = NSTEP * NCHIP

fig = plt.figure(figsize=(1.77, 1.00))   # 4.5 cm x 2.54 cm
ax_v = fig.add_axes([0.00, 0.66, 0.70, 0.30])
ax_i = fig.add_axes([0.00, 0.34, 0.70, 0.28])
ax_p = fig.add_axes([0.00, 0.02, 0.70, 0.28])
ax_s = fig.add_axes([0.76, 0.14, 0.24, 0.72])

# --- V_HCG: the wavelength staircase, steps numbered ----------------------
ax_v.step(t, volt, where='post', color=VOLT, lw=1.3)
for m in range(NSTEP):
    ax_v.text((m + 0.5) * NCHIP, m + 0.55, r'$\lambda_%d$' % (m + 1),
              fontsize=FS_SML, color=VOLT, ha='center', va='bottom')
ax_v.set_ylim(-0.6, NSTEP + 0.9)
ax_v.text(0.0, NSTEP + 0.95, r'$V_\mathrm{HCG}$: tuning', fontsize=FS_LBL,
          color=VOLT, va='bottom', ha='left')

# --- I: constant bias, the code on top -------------------------------------
ax_i.axhline(0.0, color='0.72', lw=0.6, ls=(0, (2.5, 2)))
ax_i.step(t, curr, where='post', color=CURR, lw=1.0)
ax_i.set_ylim(-0.75, 2.5)
ax_i.text(0.0, 2.45, r'$I$: bias $+$ code', fontsize=FS_LBL, color=CURR,
          va='top', ha='left')
ax_i.text(-0.3, -0.05, 'bias', fontsize=FS_SML, color='0.55', va='top',
          ha='left')
# c(t) marked over the last step, a bracket spanning one code period
x0, x1 = 3 * NCHIP, 4 * NCHIP
ax_i.plot([x0, x0, x1, x1], [1.25, 1.42, 1.42, 1.25], color=CURR, lw=0.6)
ax_i.text(0.5 * (x0 + x1), 1.5, r'$c(t)$', fontsize=FS_LBL, color=CURR,
          ha='center', va='bottom')

# --- P_out: light at lambda_m, intensity follows c(t) ----------------------
ax_p.fill_between(t, 0, curr, step='post', color=LIGHT, alpha=0.18, lw=0)
ax_p.step(t, curr, where='post', color=LIGHT, lw=0.9)
ax_p.set_ylim(-0.75, 2.5)
ax_p.text(0.0, 2.45, r'$P_\mathrm{out}$: light at $\lambda_m$ carries $c(t)$',
          fontsize=FS_LBL, color=LIGHT, va='top', ha='left')
ax_p.annotate('', xy=(T + 0.5, -0.55), xytext=(T - 6, -0.55),
              arrowprops=dict(arrowstyle='-|>', lw=0.6, color='0.45',
                              mutation_scale=5))
ax_p.text(T - 7, -0.5, r'$t$', fontsize=FS_SML, color='0.45', ha='right',
          va='center')

# step boundaries through the three traces
for a in (ax_v, ax_i, ax_p):
    for m in range(1, NSTEP):
        a.axvline(m * NCHIP, color='0.82', lw=0.4, ls=(0, (1, 1.5)), zorder=0)
    a.set_xlim(-0.5, T + 0.5)
    a.axis('off')

# --- laser spectrum: one narrow line per step ------------------------------
lam = np.linspace(-0.6, NSTEP - 0.4, 1200)
w = 0.045
for m in range(NSTEP):
    line = 1.0 / (1.0 + ((lam - m) / w) ** 2)
    ax_s.plot(lam, line, color=VOLT, lw=0.8)
    ax_s.text(m, -0.12, r'$\lambda_%d$' % (m + 1), fontsize=4.6, color=VOLT,
              ha='center', va='top')
ax_s.annotate('', xy=(NSTEP - 0.35, 0.0), xytext=(-0.6, 0.0),
              arrowprops=dict(arrowstyle='-|>', lw=0.6, color='0.45',
                              mutation_scale=5))
ax_s.text(NSTEP - 0.3, 0.02, r'$\lambda$', fontsize=FS_SML, color='0.45',
          ha='left', va='bottom')
ax_s.set_xlim(-0.7, NSTEP - 0.1)
ax_s.set_ylim(-0.42, 1.55)
ax_s.text(-0.6, 1.5, 'spectrum', fontsize=FS_LBL, color=VOLT, va='top',
          ha='left')
ax_s.axis('off')

os.makedirs('figs/fig1_panels', exist_ok=True)
fig.savefig('figs/fig1_panels/inset_drive.pdf', bbox_inches='tight',
            pad_inches=0.01, transparent=True)
fig.savefig('figs/fig1_panels/inset_drive.png', dpi=400, bbox_inches='tight',
            pad_inches=0.01)
fig.savefig('figs/fig1_panels/inset_drive.svg', bbox_inches='tight',
            pad_inches=0.01, transparent=True)
plt.close(fig)
print('saved figs/fig1_panels/inset_drive.pdf, .png, .svg')
