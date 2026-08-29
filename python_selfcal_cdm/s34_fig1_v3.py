"""s34_fig1_v3.py - alternative layout for Fig. 1, for the author to judge.

Not wired into the paper. main.tex still uses the TikZ composition. This is
the layout the reviewer proposed on 29.08, built so that the choice can be
made on a rendered figure rather than on a description.

What a first-time reader has to get in ten seconds, in this order:
  1. one beam carries a code and steps through wavelength;
  2. every grating answers in the SAME band and differs only by delay;
  3. correlation turns one mixed trace into a delay-wavelength map;
  4. the product is one number, the fitted Bragg wavelength, read off a CUT
     through that map, not the map itself.

Layout: a full-width strip that answers "where", then three panels that
answer "what the receiver does". Two arrows, one word each. Three gratings
in three saturated colors, everything else gray. Dropped from the current
figure: circled step numbers, the raw record as a panel of its own, the
staircase drive plot, the overlapped-spectra panel without axes, and the
note 'color = position, not wavelength', which stops being needed once
color means grating and nothing else.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec
import common as C
import figstyle as FS

FS.apply(7.0)

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
COLS = [FS.BLUE, FS.ORANGE, FS.GREEN]          # one color per grating
GREY = '0.55'

# ---------------------------------------------------------------------------
# physics: three gratings, one band, separated by round-trip delay only
# ---------------------------------------------------------------------------
N = 31                                          # short code, so chips are visible
code = C._mls01(5).astype(float)
DELAYS = [4, 13, 23]                            # chips
AMPS = [1.00, 0.62, 0.38]                       # thinning along the fiber
DETUNE = [-90.0, 10.0, 120.0]                   # pm, all inside one band

M = 41
lam = np.linspace(-260, 260, M)                 # pm around the band centre
nu = lam * C.GHZ_PER_PM

record = np.zeros(N)
for d, a, dl in zip(DELAYS, AMPS, DETUNE):
    record += a * np.roll(code, d) * 0.9

rep = 2.0 * code - 1.0                          # bipolar replica, sum +1
sweep = np.zeros((M, N))
for m in range(M):
    row = np.zeros(N)
    for d, a, dl in zip(DELAYS, AMPS, DETUNE):
        row += a * C.fbg_gauss(nu[m], dl * C.GHZ_PER_PM, F) * np.roll(code, d)
    sweep[m] = row
corr = np.array([np.real(np.fft.ifft(np.fft.fft(r) * np.fft.fft(rep).conj()))
                 for r in sweep])
corr /= corr.max()

cut = corr[:, DELAYS[1]]
fit_pm = C.gauss_fit_peak(nu, cut) * PM

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(7.1, 3.35))
gs = GridSpec(2, 3, height_ratios=[0.46, 1.0], hspace=0.30, wspace=0.34,
              left=0.065, right=0.985, top=0.99, bottom=0.115)

# --- strip: where the gratings are ----------------------------------------
top = fig.add_subplot(gs[0, :])
top.set_xlim(0, 10); top.set_ylim(-0.55, 1.15); top.axis('off')
top.plot([1.45, 9.7], [0.45, 0.45], color='0.35', lw=1.3)
top.add_patch(Rectangle((0.05, 0.24), 0.92, 0.42, fc='white', ec='0.35', lw=0.8))
top.text(0.51, 0.45, 'swept\nVCSEL', ha='center', va='center', fontsize=6.4)
top.plot([1.16], [0.45], marker='o', ms=7, mfc='white', mec='0.35', mew=0.8)
top.text(1.16, 0.72, 'circ.', ha='center', fontsize=5.6, color='0.45')
top.annotate('', xy=(1.16, -0.14), xytext=(1.16, 0.30),
             arrowprops=dict(arrowstyle='-|>', color='0.35', lw=0.9))
top.text(1.16, -0.30, 'PD', ha='center', va='center', fontsize=6.4)
for i, (z, col) in enumerate(zip([3.6, 6.1, 8.7], COLS)):
    for off in (-0.05, 0.0, 0.05):
        top.plot([z + off, z + off], [0.30, 0.60], color=col, lw=1.5)
    top.text(z, 0.74, r'$z_%d$' % (i + 1), ha='center', fontsize=7.0, color=col)
    top.text(z, 0.10, r'$\tau_%d$' % (i + 1), ha='center', fontsize=7.0, color=col)
top.text(9.7, -0.28, r'$\tau_k = 2n_g z_k/c$', ha='right', fontsize=7.2,
         color='0.35')
top.text(5.6, -0.30, 'the round trip is the only thing that tells them apart',
         ha='center', fontsize=6.2, color='0.45', style='italic')

# --- (a) one code period at one wavelength step ---------------------------
a = fig.add_subplot(gs[1, 0])
t = np.arange(N)
a.step(t, 3.35 + 0.72 * code, where='post', color='0.25', lw=1.0)
a.text(0.4, 4.22, 'transmitted code', fontsize=6.0, color='0.25')
for i, (d, amp, col) in enumerate(zip(DELAYS, AMPS, COLS)):
    y0 = 2.20 - 0.72 * i
    a.step(t, y0 + 0.46 * amp * np.roll(code, d), where='post', color=col, lw=0.9)
a.step(t, -0.62 + 0.42 * record / record.max(), where='post', color=GREY, lw=0.9)
a.text(0.4, -0.05, 'their sum, what the PD sees', fontsize=6.0, color=GREY)
a.annotate('', xy=(DELAYS[1], 2.98), xytext=(0, 2.98),
           arrowprops=dict(arrowstyle='<->', color=COLS[1], lw=0.7))
a.text(0.5, 3.05, r'delayed by $\tau_2$', ha='left', va='bottom',
       fontsize=6.0, color=COLS[1])
a.set_xlim(0, N); a.set_ylim(-0.80, 4.55)
a.set_yticks([]); a.set_xlabel('time [chips]')
a.set_title('(a) one code period, one wavelength step', fontsize=7.4, loc='left')
for sp in ('top', 'right', 'left'):
    a.spines[sp].set_visible(False)

# --- (b) the map ----------------------------------------------------------
b = fig.add_subplot(gs[1, 1])
b.imshow(corr, aspect='auto', origin='lower', cmap='Greys',
         extent=[0, N, lam[0], lam[-1]], vmin=0, vmax=1)
for d, dl, col in zip(DELAYS, DETUNE, COLS):
    b.plot([d + 0.5], [dl], marker='v', ms=4.5, color=col, clip_on=False)
b.axvline(DELAYS[1] + 0.5, color=COLS[1], ls=(0, (2, 1.6)), lw=0.9)
b.set_xlim(0, N); b.set_ylim(lam[0], lam[-1])
b.set_xlabel(r'delay $\tau$ [chips]')
b.set_ylabel(r'wavelength offset [pm]')
b.set_title('(b) after correlating every row', fontsize=7.4, loc='left')

# --- (c) the cut, which is the measurement --------------------------------
c = fig.add_subplot(gs[1, 2])
c.plot(lam, cut / cut.max(), 'o', ms=2.6, color=COLS[1], mfc='white', mew=0.7,
       label='samples')
fine = np.linspace(lam[0], lam[-1], 400)
c.plot(fine, C.fbg_gauss(fine * C.GHZ_PER_PM, fit_pm * C.GHZ_PER_PM, F),
       color=COLS[1], lw=1.1, label='fit')
c.axvline(fit_pm, color='0.3', ls=(0, (2, 1.6)), lw=0.8)
c.annotate(r'$\hat\lambda_{B,2}$', xy=(fit_pm, 1.06), xytext=(fit_pm + 55, 1.14),
           fontsize=7.2, color='0.25',
           arrowprops=dict(arrowstyle='-', color='0.3', lw=0.6))
c.annotate('', xy=(fit_pm + 205, 0.40), xytext=(fit_pm + 145, 0.40),
           arrowprops=dict(arrowstyle='-|>', color='0.45', lw=0.8))
c.text(fit_pm + 175, 0.45, r'$\Delta T,\ \varepsilon$', fontsize=6.4,
       color='0.45', ha='center', va='bottom')
c.set_xlim(lam[0], lam[-1]); c.set_ylim(0, 1.30)
c.set_xlabel(r'wavelength offset [pm]'); c.set_ylabel('reflectivity')
c.set_title('(c) the cut at $\\tau_2$ is the measurement', fontsize=7.4,
            loc='left')
c.legend(fontsize=5.8, loc='upper left', frameon=False, handlelength=1.2)

fig.savefig('figs/fig_s34_concept_v3.pdf', bbox_inches='tight', pad_inches=0.02)
fig.savefig('figs/fig_s34_concept_v3.png', dpi=300, bbox_inches='tight',
            pad_inches=0.02)
plt.close(fig)

print('alternative Fig. 1 layout, for comparison against the TikZ version')
print('  gratings at chips %s, detunings %s pm' % (DELAYS, DETUNE))
print('  cut at tau_2 fits to %.1f pm against a true %.1f pm'
      % (fit_pm, DETUNE[1]))
print('saved figs/fig_s34_concept_v3.png and .pdf')
