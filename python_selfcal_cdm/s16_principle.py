"""s16_principle.py - the two explanatory figures of the paper.

Figure 1, principle: the optical layout on the left and, on the right, what the
correlator actually returns. The despread output is a two-dimensional map, delay
bin against wavelength: the code marks WHICH wavelength was launched, the delay
marks WHICH grating sent the light back. A horizontal cut through that map is one
grating's reflection spectrum, and its peak is the measurement.

Figure 2, mechanisms: the three array effects that this study quantifies, each
drawn from the same model used for the results, not sketched by hand.
  (a) multiple reflections: a third-order path returns at tau_a - tau_b + tau_c,
      which for uniform spacing is always another occupied bin and for randomised
      spacing usually is not;
  (b) code leakage: the autocorrelation side lobe adds a scaled copy of every
      other grating's spectrum as a broad background under the wanted line.

Spectral shadowing, the third array term, has its own figure in s19 because it
comes with a correction.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS
FS.apply()

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
NCH = 127
MSEQ = 1.0 - 2.0 * C._mls01(7)
ACORR = C.periodic_xcorr(MSEQ, MSEQ)


# ===========================================================================
# FIGURE 1: principle
# ===========================================================================
fig = plt.figure(figsize=(7.1, 2.65))
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.28)

# --- (a) optical layout ----------------------------------------------------
ax = fig.add_subplot(gs[0, 0]); ax.axis('off')
ax.set_xlim(0, 10.4); ax.set_ylim(0, 6.2)


def box(x, y, w, h, t, fc='#eaf1fb', fs=7.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle='round,pad=0.05,rounding_size=0.12',
                                fc=fc, ec='#1f4e79', lw=1.2))
    ax.text(x + w / 2, y + h / 2, t, ha='center', va='center', fontsize=fs)


def arrow(x1, y1, x2, y2, col='#333', lw=1.2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>',
                                 mutation_scale=10, lw=lw, color=col))


box(0.15, 4.20, 2.55, 1.35, 'swept VCSEL,\ncode-modulated', fs=7.4)
# code waveform inside the source box area
tt = np.linspace(0.35, 2.35, 300)
bits = MSEQ[:16]
cw = np.repeat(bits, len(tt) // 16 + 1)[:len(tt)]
ax.plot(tt, 3.55 + 0.18 * cw, color='#1f4e79', lw=0.9)
ax.text(1.35, 3.15, 'code c(t)', ha='center', fontsize=7, color='#1f4e79')

ax.add_patch(Circle((3.5, 4.85), 0.42, fc='#fdf3e7', ec='#1f4e79', lw=1.2))
ax.text(3.5, 4.85, 'circ', ha='center', va='center', fontsize=7)
arrow(2.72, 4.85, 3.05, 4.85)

# fiber with gratings at irregular positions
ax.plot([3.95, 10.2], [4.85, 4.85], color='#555', lw=1.6)
zpos = [4.7, 5.9, 7.6, 8.5, 9.7]
for i, z in enumerate(zpos):
    ax.add_patch(Rectangle((z - 0.16, 4.62), 0.32, 0.46, fc='#c0392b',
                           ec='#7b241c', hatch='///', lw=0.8))
    ax.text(z, 5.30, r'$z_%d$' % (i + 1), ha='center', fontsize=7)
ax.annotate('', xy=(4.7, 4.35), xytext=(5.9, 4.35),
            arrowprops=dict(arrowstyle='<->', lw=0.8, color='#555'))
ax.annotate('', xy=(7.6, 4.35), xytext=(8.5, 4.35),
            arrowprops=dict(arrowstyle='<->', lw=0.8, color='#555'))
ax.text(7.2, 3.95, 'randomized spacing (Section IV-B)', ha='center', fontsize=7,
        color='#555')
ax.text(7.1, 5.75, 'FBG array: same nominal wavelength,\nseparated only by delay',
        ha='center', fontsize=7.4)

# return path
arrow(3.5, 4.43, 3.5, 3.0)
box(2.75, 2.05, 1.5, 0.9, 'InGaAs\nAPD')
arrow(4.25, 2.5, 5.05, 2.5)
box(5.00, 2.05, 2.45, 0.9, 'equivalent-time\nsampling', fs=7.2)
arrow(7.45, 2.5, 8.15, 2.5)
box(8.15, 2.05, 2.2, 0.9, 'matched-filter\ncorrelator', fs=7.2)
arrow(9.25, 2.05, 9.25, 1.35)
box(6.6, 0.35, 3.7, 1.0, 'delay-wavelength map  ->  peak fit  ->  $\\lambda_B$',
    fc='#eafaf0')

ax.text(0.1, 1.3, 'one code on the fiber\nat a time: the sweep\nvisits wavelengths\n'
                  'in sequence, so gratings\nare separated by the\nAUTOCORRELATION\n'
                  'side lobe', fontsize=7.2, va='center', color='#1f4e79')
ax.set_title('(a) Coded interrogation of a grating array', fontsize=9)

# --- (b) despread map ------------------------------------------------------
ax2 = fig.add_subplot(gs[0, 1])
M = 96
nu = np.linspace(-2.6 * F, 2.6 * F, M)
rng = np.random.default_rng(2)
K = 5
bins = np.array([11, 27, 46, 79, 103])
nub = np.array([-18.0, 9.0, -4.0, 21.0, -11.0])
R = 0.05
shapes = np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
tcum = np.ones((K, M))
for k in range(1, K):
    tcum[k] = tcum[k - 1] * (1.0 - R * shapes[k - 1]) ** 2
prim = R * shapes * tcum
grid = np.zeros((NCH, M))
grid[bins] = prim
W = ACORR[(np.arange(NCH)[:, None] - bins[None, :]) % NCH]
grid = grid + W @ prim
im = ax2.imshow(grid.T, aspect='auto', origin='lower', cmap='viridis',
                extent=[0, NCH, nu[0] * PM / 1000.0, nu[-1] * PM / 1000.0])
for k in range(K):
    ax2.plot(bins[k], nub[k] * PM / 1000.0, 'o', mfc='none', mec='w', ms=9, mew=1.2)
    ax2.text(bins[k] + 2, nub[k] * PM / 1000.0 + 0.06, 'FBG %d' % (k + 1),
             color='w', fontsize=7)
ax2.set_xlabel('delay bin  [chip]  ->  grating position')
ax2.set_ylabel('wavelength offset [nm]')
ax2.set_title('(b) What the correlator returns', fontsize=9)
cb = fig.colorbar(im, ax=ax2, pad=0.02); cb.set_label('despread reflectance', fontsize=7.5)
cb.ax.tick_params(labelsize=7)
ax2.text(0.98, 0.04, 'a horizontal cut is one grating spectrum',
         transform=ax2.transAxes, ha='right', fontsize=7, color='w')
fig.savefig('figs/fig_s16_principle.png', dpi=150, bbox_inches='tight')

# ===========================================================================
# FIGURE 2: error mechanisms
# ===========================================================================
# This figure is included at 0.72 text width in the paper. Draw it at that
# physical width so the 7 pt labels remain 7 pt after LaTeX placement.
fig2, ax = plt.subplots(1, 2, figsize=(5.1, 2.30))

# --- (b) ghost delays ------------------------------------------------------
axb = ax[0]
axb.set_xlim(0, 10); axb.set_ylim(-1.6, 10); axb.axis('off')
axb.plot([0.6, 9.4], [8.6, 8.6], color='#555', lw=1.4)
gz = [2.0, 4.4, 7.4]
for z, nm in zip(gz, ['b', 'a', 'c']):
    axb.add_patch(Rectangle((z - 0.16, 8.35), 0.32, 0.5, fc='#c0392b',
                            ec='#7b241c', hatch='///', lw=0.7))
    axb.text(z, 9.15, nm, ha='center', fontsize=8)
axb.annotate('', xy=(4.4, 7.9), xytext=(0.7, 7.9),
             arrowprops=dict(arrowstyle='-|>', color='#28b463', lw=1.1))
axb.annotate('', xy=(2.0, 7.35), xytext=(4.4, 7.35),
             arrowprops=dict(arrowstyle='-|>', color='#28b463', lw=1.1))
axb.annotate('', xy=(7.4, 6.8), xytext=(2.0, 6.8),
             arrowprops=dict(arrowstyle='-|>', color='#28b463', lw=1.1))
axb.annotate('', xy=(0.7, 6.25), xytext=(7.4, 6.25),
             arrowprops=dict(arrowstyle='-|>', color='#28b463', lw=1.1))
axb.text(5.0, 5.55, r'ghost: $\tau_a-\tau_b+\tau_c$, power $\propto R^3$',
         ha='center', fontsize=7.0)

Kg = 10
bu = 3 + 7 * np.arange(Kg)
rng2 = np.random.default_rng(11)
br = np.sort(rng2.choice(np.arange(1, 90), size=Kg, replace=False))


def ghost_bins(b):
    out = []
    for a in range(len(b)):
        for bb in range(len(b)):
            for c in range(len(b)):
                if bb < a and bb < c:
                    gb = b[a] - b[bb] + b[c]
                    if b.min() <= gb <= b.max():
                        out.append(gb)
    return np.array(out)


for row, (b, gbs, ttl, col) in enumerate(
        [(bu, ghost_bins(bu), 'uniform spacing', '#c0392b'),
         (br, ghost_bins(br), 'randomized spacing', '#2980b9')]):
    y0 = 3.3 - 2.6 * row
    axb.plot([0.6, 9.4], [y0, y0], color='0.8', lw=0.8)
    xs = 0.6 + 8.8 * (b - 0) / 90.0
    axb.plot(xs, [y0] * len(b), 'v', color=col, ms=6, label='gratings' if row == 0 else None)
    hit = np.array([g in b for g in gbs])
    xg = 0.6 + 8.8 * gbs / 90.0
    axb.plot(xg[hit], [y0 + 0.55] * hit.sum(), 'x', color='#c0392b', ms=4.5, mew=1.0)
    axb.plot(xg[~hit], [y0 + 0.55] * (~hit).sum(), '.', color='0.6', ms=3)
    axb.text(0.6, y0 - 0.75, '%s: %d%% collide'
             % (ttl, round(100 * hit.mean())), fontsize=7.0, color=col)
axb.text(5.0, -1.15, 'delay bin', ha='center', fontsize=7.2, color='0.35')
axb.set_title('(a) Spacing decides ghost collisions', fontsize=8.1)

# --- (c) code leakage ------------------------------------------------------
Mc = 96
nuc = np.linspace(-2.6 * F, 2.6 * F, Mc)
rng3 = np.random.default_rng(7)
Kc = 32
nubs = rng3.uniform(-25, 25, Kc)
Ac = 0.05 * np.exp(-0.5 * ((nuc[None, :] - nubs[:, None]) / SIG) ** 2)
wanted = Ac[0]
leak = (-1.0 / NCH) * (Ac[1:].sum(axis=0))
ax[1].plot(nuc * PM / 1000.0, wanted, color='#2980b9', lw=1.4, label='wanted grating')
ax[1].plot(nuc * PM / 1000.0, np.abs(leak), color='0.45', lw=1.2,
           label=r'|leakage|, side lobe $1/N$')
ax[1].plot(nuc * PM / 1000.0, np.abs((-17.0 / NCH) * (Ac[1:].sum(axis=0))),
           color='#c0392b', lw=1.2, ls='--', label=r'|leakage|, Gold $17/N$')
ax[1].set_yscale('log'); ax[1].set_ylim(1e-5, 0.2)
ax[1].set_xlabel('wavelength offset [nm]'); ax[1].set_ylabel('despread reflectance')
ax[1].set_title('(b) Side-lobe leakage, K = %d' % Kc, fontsize=8.1)
handles, labels = ax[1].get_legend_handles_labels()
ax[1].grid(True, which='both', alpha=0.2)

fig2.subplots_adjust(left=0.055, right=0.99, top=0.88, bottom=0.36,
                     wspace=0.30)
fig2.legend(handles, labels, fontsize=6.3, loc='lower center',
            bbox_to_anchor=(0.75, 0.015), frameon=False, ncol=1,
            handlelength=2.3, labelspacing=0.16)
fig2.savefig('figs/fig_s16_mechanisms.png', dpi=150, bbox_inches='tight')
fig2.savefig('figs/fig_s16_mechanisms.pdf', bbox_inches='tight')

print('ghost bins, uniform: %d in span, %.0f%% on a grating bin'
      % (len(ghost_bins(bu)), 100 * np.mean([g in bu for g in ghost_bins(bu)])))
print('ghost bins, random : %d in span, %.0f%% on a grating bin'
      % (len(ghost_bins(br)), 100 * np.mean([g in br for g in ghost_bins(br)])))
print('saved figs/fig_s16_principle.png, figs/fig_s16_mechanisms.png')
