"""s24_concept.py - data panels for the concept figure (Fig. 1), v9.

The figure itself is composed in TikZ (Draft_v2/fig1_concept.tex); this
script renders only the DATA panels, each as a tight standalone PDF+PNG,
with STIX (Times-like) typography so they match the IEEEtran page:

  panel_record  the raw photodiode record: M sweep rows x N chips of summed
                delayed codes, with the wavelength staircase strip at left
                (row index IS wavelength). No grating visible by eye.
  panel_map     the same record after correlating every row with c(t):
                K clean columns, one per grating (full model: shadowing
                + code leakage).
  panel_spec    one column read out: the sampled spectrum, dots per step,
                Gaussian-fitted lambda_B.
  inset_drive   staircase drive, one code period per step.
  inset_echoes  the K delayed, weighted copies of the code.
  inset_overlap the K overlapping reflection spectra.
  inset_template  one period of the template c(t).

Model identical to the article's Section II: m-sequence N=127, K=4 at
R=5%, shadowing via cumulative transmission, leakage via the periodic
autocorrelation. Nothing is drawn by hand: every panel is computed.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import common as C

# --- Times-like typography to match IEEEtran --------------------------------
plt.rcParams.update({
    'font.size': 6.0,
    'font.family': 'serif',
    'font.serif': ['STIXGeneral', 'Times New Roman', 'DejaVu Serif'],
    'mathtext.fontset': 'stix',
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'xtick.major.size': 1.8, 'ytick.major.size': 1.8,
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.dpi': 200, 'savefig.dpi': 600,
})

BLUE = '#0072B2'
ECOL = ['#D55E00', '#E69F00', '#009E73', '#CC79A7']

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
NCH = 127
_MS = 1.0 - 2.0 * C._mls01(7)
_MS = np.roll(_MS, -int(np.argmax(_MS == 1.0)))
ACORR = C.periodic_xcorr(_MS, _MS)

# --- model ------------------------------------------------------------------
M = 64
nu = np.linspace(-2.3 * F, 2.3 * F, M)
lam = nu * PM / 1000.0
K = 4
bins = np.array([16, 47, 76, 105])
nub = np.array([-16.0, 9.0, -5.0, 19.0])
R = 0.05
shapes = np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
tcum = np.ones((K, M))
for k in range(1, K):
    tcum[k] = tcum[k - 1] * (1.0 - R * shapes[k - 1]) ** 2
prim = R * shapes * tcum
Wl = ACORR[(np.arange(NCH)[:, None] - bins[None, :]) % NCH]
corr = (Wl @ prim).T                     # (M, NCH) full model
HL = 1
JM = np.argmin(np.abs(nu - 0.0))
spec = corr[:, bins[HL]]
mu_fit = C.gauss_fit_peak(nu, spec)
mu_row = np.interp(mu_fit, nu, np.arange(M))

OS = 2
code01 = np.repeat(0.5 * (1 + _MS), OS)
rng = np.random.default_rng(7)
rec = np.zeros((M, NCH * OS))
for k in range(K):
    rec += prim[k][:, None] * np.roll(code01, bins[k] * OS)[None, :]
rec += rng.normal(0, 0.004, rec.shape)

OUT = 'figs/fig1_panels/'
import os
os.makedirs(OUT, exist_ok=True)


def save(fig, name):
    fig.savefig(OUT + name + '.pdf', bbox_inches='tight', pad_inches=0.01)
    fig.savefig(OUT + name + '.png', bbox_inches='tight', pad_inches=0.01)
    plt.close(fig)
    print('saved', OUT + name)


# --- panel_record: staircase strip + raw record -----------------------------
fig = plt.figure(figsize=(2.30, 1.95))
axl = fig.add_axes([0.135, 0.165, 0.135, 0.815])
axb = fig.add_axes([0.300, 0.165, 0.690, 0.815])
axl.plot(lam, np.arange(M), color=BLUE, lw=0.7, drawstyle='steps-mid')
axl.plot([lam[JM]], [JM], marker='o', ms=2.0, color='0.15', zorder=5)
axl.set_ylim(-0.5, M - 0.5); axl.set_xlim(-0.55, 0.55)
axl.set_xticks([]); axl.set_yticks([0, M - 1])
axl.set_yticklabels(['1', '$M$'])
axl.tick_params(labelsize=5.5, pad=1.2)
axl.set_xlabel(r'$\lambda_m$', fontsize=6, labelpad=1)
axl.set_ylabel('sweep step $m$', fontsize=6, labelpad=0)
axl.spines['bottom'].set_visible(False)
axb.imshow(rec, aspect='auto', origin='lower', cmap='Greys',
           extent=[0, NCH, -0.5, M - 0.5], interpolation='nearest',
           vmax=rec.max())
axb.axhline(JM, color=BLUE, lw=0.5, ls=(0, (3, 1.5)))
axb.text(4, JM + 2.0, r'step $\lambda_m$', fontsize=5.5, color=BLUE,
         bbox=dict(fc='white', ec='none', alpha=0.75, pad=0.5))
axb.set_xlim(0, NCH); axb.set_ylim(-0.5, M - 0.5)
axb.set_xticks([0, NCH]); axb.set_xticklabels(['0', '$N$'])
axb.set_yticks([])
axb.tick_params(labelsize=5.5, pad=1.2)
axb.set_xlabel('time in code period (chips)', fontsize=6, labelpad=1)
save(fig, 'panel_record')

# --- panel_map: correlated record -------------------------------------------
fig = plt.figure(figsize=(2.10, 1.95))
ax = fig.add_axes([0.06, 0.165, 0.90, 0.815])
ax.imshow(corr, aspect='auto', origin='lower', cmap='Greys',
          extent=[0, NCH, -0.5, M - 0.5], vmax=corr.max() * 0.9)
for b, col in zip(bins, ECOL):
    ax.plot(b, M + 0.6, marker='v', ms=2.6, color=col, clip_on=False)
ax.axvspan(bins[HL] - 3.5, bins[HL] + 3.5, color=ECOL[HL], alpha=0.16, lw=0)
ax.axhline(JM, color=BLUE, lw=0.5, ls=(0, (3, 1.5)))
ax.text(123, JM + 2.0, 'same row', fontsize=5.5, color=BLUE, ha='right')
ax.text(bins[HL] + 5.5, 1.5, 'column $=$ grating', fontsize=5.5,
        color='0.15', ha='left', va='bottom')
ax.set_xlim(0, NCH); ax.set_ylim(-0.5, M - 0.5)
ax.set_xticks(list(bins))
ax.set_xticklabels([r'$\tau_1$', r'$\tau_2$', r'$\tau_3$', r'$\tau_4$'],
                   fontsize=6)
for t, col in zip(ax.get_xticklabels(), ECOL):
    t.set_color(col)
ax.set_yticks([])
ax.tick_params(labelsize=5.5, pad=1.2)
ax.set_xlabel(r'delay $\tau$ (chips)', fontsize=6, labelpad=1)
save(fig, 'panel_map')

# --- panel_spec: one column read out ----------------------------------------
fig = plt.figure(figsize=(1.45, 1.95))
ax = fig.add_axes([0.10, 0.165, 0.86, 0.815])
sn = spec / spec.max()
ax.plot(sn, np.arange(M), color=ECOL[HL], lw=0.9)
ax.plot(sn[::2], np.arange(M)[::2], ls='none', marker='o', ms=1.5,
        color=ECOL[HL])
ax.axhline(mu_row, color='0.35', ls=(0, (3, 2)), lw=0.6)
ax.text(1.00, mu_row + 1.8, r'$\lambda_B$', fontsize=7, color='0.15',
        ha='right')
ax.text(1.02, 4.0, r'$\Delta\lambda_B \propto \Delta T,\ \varepsilon$',
        fontsize=6, color='0.15', ha='right')
ax.set_xlim(-0.04, 1.06); ax.set_ylim(-0.5, M - 0.5)
ax.set_xticks([0, 1]); ax.set_yticks([])
ax.tick_params(labelsize=5.5, pad=1.2)
ax.set_xlabel('reflectivity', fontsize=6, labelpad=1)
save(fig, 'panel_spec')

# --- insets (frameless content; TikZ draws the frames) ----------------------
OSR = 6
NSTEP, CHPS = 4, 14
stair = []
for i in range(NSTEP):
    stair.append(i * 1.0 + 0.42 * 0.5 * (1 + _MS[i * CHPS:(i + 1) * CHPS]))
stair = np.repeat(np.concatenate(stair), OSR)
ts = np.arange(len(stair)) / (OSR * CHPS)

fig = plt.figure(figsize=(0.95, 0.72))
ax = fig.add_axes([0.02, 0.05, 0.96, 0.92])
ax.plot(ts, stair, color=BLUE, lw=0.6)
ax.set_xlim(-0.05, ts[-1] + 0.05); ax.set_ylim(-0.15, NSTEP + 0.1)
ax.axis('off')
save(fig, 'inset_drive')

taus = [3.0, 7.5, 12.0, 16.5]
drv2 = np.repeat(0.5 * (1 + _MS[:40]), OSR)
tt2 = np.arange(len(drv2)) / OSR
W2 = tt2 < 34
fig = plt.figure(figsize=(1.15, 0.55))
ax = fig.add_axes([0.02, 0.05, 0.96, 0.92])
for tcau, a, col, off in zip(taus, [1.0, 0.80, 0.62, 0.48], ECOL,
                             (3.6, 2.4, 1.2, 0.0)):
    e = np.zeros_like(drv2); sh = int(tcau * OSR)
    e[sh:] = a * drv2[:len(drv2) - sh]
    ax.plot(tt2[W2], off + e[W2], color=col, lw=0.5)
ax.set_xlim(0, 34); ax.set_ylim(-0.3, 5.0)
ax.axis('off')
save(fig, 'inset_echoes')

fig = plt.figure(figsize=(0.95, 0.62))
ax = fig.add_axes([0.02, 0.06, 0.96, 0.90])
lamg = np.linspace(-60, 60, 200)
for c0, col in zip(nub, ECOL):
    g = np.exp(-0.5 * ((lamg - c0) / SIG) ** 2)
    ax.plot(lamg, g, color=col, lw=0.55)
    ax.fill_between(lamg, g, color=col, alpha=0.10, lw=0)
ax.set_xlim(-60, 60); ax.set_ylim(0, 1.08)
ax.axis('off')
save(fig, 'inset_overlap')

fig = plt.figure(figsize=(0.85, 0.34))
ax = fig.add_axes([0.02, 0.10, 0.96, 0.80])
ax.plot(np.arange(len(code01)) / OS, code01, color='0.3', lw=0.55)
ax.set_xlim(0, 22); ax.set_ylim(-0.25, 1.25)
ax.axis('off')
save(fig, 'inset_template')

print('fitted lambda_B: %.1f pm (true %.1f pm)' % (mu_fit * PM, nub[HL] * PM))
