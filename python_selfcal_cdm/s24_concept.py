"""s24_concept.py - data panels for the concept figure (Fig. 1), v11.

Redesign after a three-critic review (graphics editor, first-glance
physicist, information designer). Two visual levels only:

  level 1  a heavy hardware spine (drawn in TikZ), with two computed
           micro-insets: the coded staircase drive, and the delayed echoes
           with their black SUM and a dimensioned tau_2, plus the K
           overlapping spectra crowded around one lambda_B.
  level 2  ONE hero: the delay-wavelength map as a ridgeline plot, each
           grating's ridge in its own colour (the JLT Fig. 1(e) lineage),
           flanked by a stamp-sized raw record ("no structure") on the left
           and the extracted spectrum on the right, sharing the lambda axis.

Model as in the paper: m-sequence N=127, K=4 at R=5%, shadowing via
cumulative transmission, leakage via the periodic autocorrelation. M=24
sweep steps for a readable ridge density.

Panels written to figs/fig1_panels/:
  panel_ridge   the hero ridgeline map (delay x, wavelength y)
  panel_spec    extracted spectrum, vertical lambda axis matched to ridge
  panel_stamp   small raw record, the "before" swatch
  inset_drive   staircase drive, one code period per step (square steps)
  inset_echosum four delayed, clearly weighted echoes + their black sum,
                with a tau_2 dimension against the transmit marker
  inset_overlap K near-identical spectra crowded around one lambda_B
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings, os; warnings.filterwarnings('ignore')
import common as C

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
ECOL = ['#C74E0A', '#E8A200', '#009E73', '#CC79A7']   # tau1 shifted redward,
                                                      # tau2 toward gold
PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
NCH = 127
_MS = 1.0 - 2.0 * C._mls01(7)
_MS = np.roll(_MS, -int(np.argmax(_MS == 1.0)))
ACORR = C.periodic_xcorr(_MS, _MS)

# --- model ------------------------------------------------------------------
M = 24
nu = np.linspace(-2.05 * F, 2.05 * F, M)
lam = nu * PM / 1000.0                       # nm offsets
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
corr = (Wl @ prim).T                         # (M, NCH)
# render the two-chip peak width: triangular smear over delay
kern = np.array([0.45, 1.0, 0.45]); kern /= kern.max()
corr_d = np.apply_along_axis(lambda r: np.convolve(r, kern, 'same'), 1, corr)
corr_d = corr_d / (R)                        # ~unit peak scale
HL = 1
spec = corr[:, bins[HL]] / corr[:, bins[HL]].max()
A_f, mu_f, sig_f, base_f = C.gauss_fit_full(nu, corr[:, bins[HL]])

OS = 2
code01 = np.repeat(0.5 * (1 + _MS), OS)
rng = np.random.default_rng(7)
rec = np.zeros((M, NCH * OS))
for k in range(K):
    rec += prim[k][:, None] * np.roll(code01, bins[k] * OS)[None, :]
rec += rng.normal(0, 0.004, rec.shape)

OUT = 'figs/fig1_panels/'
os.makedirs(OUT, exist_ok=True)


def save(fig, name):
    fig.savefig(OUT + name + '.pdf', bbox_inches='tight', pad_inches=0.01)
    fig.savefig(OUT + name + '.png', bbox_inches='tight', pad_inches=0.01)
    plt.close(fig)
    print('saved', OUT + name)


# --- shared vertical lambda frame for ridge + spec --------------------------
GAP = lam[1] - lam[0]                        # nm per step
AMP = 4.6 * GAP                              # ridge amplitude
YLO = lam[0] - 1.2 * GAP
YHI = lam[-1] + AMP + 1.2 * GAP

# =============================== hero: ridgeline ============================
fig = plt.figure(figsize=(3.80, 2.45))
ax = fig.add_axes([0.085, 0.125, 0.895, 0.860])
x = np.arange(NCH)
for m in range(M - 1, -1, -1):               # back to front
    y = lam[m] + AMP * corr_d[m]
    ax.fill_between(x, lam[m], y, color='white', lw=0, zorder=2 * m)
    ax.plot(x, y, color='0.68', lw=0.4, zorder=2 * m + 1)
    for k in range(K):
        if corr_d[m, bins[k]] < 0.055:       # colour only where a peak is
            continue
        s = slice(max(bins[k] - 3, 0), bins[k] + 4)
        ax.fill_between(x[s], lam[m], y[s], color=ECOL[k], alpha=0.55,
                        lw=0, zorder=2 * m + 1)
        ax.plot(x[s], y[s], color=ECOL[k], lw=0.8, zorder=2 * m + 1)
# envelope of the highlighted ridge
env = lam + AMP * corr_d[:, bins[HL]]
ax.plot(np.full(M, bins[HL]), env, color='0.15', lw=0.6,
        ls=(0, (2.5, 1.6)), zorder=100)
ax.annotate('one ridge $=$ one grating', xy=(bins[3] + 2, env.max()),
            xytext=(60, YHI - 0.35 * GAP), fontsize=6.2, color='0.15',
            ha='left', va='top', zorder=101)
ax.set_xlim(-2, NCH + 1)
ax.set_ylim(YLO, YHI)
ax.set_xticks(list(bins))
ax.set_xticklabels([r'$\tau_1$', r'$\tau_2$', r'$\tau_3$', r'$\tau_4$'],
                   fontsize=6.4)
for t, col in zip(ax.get_xticklabels(), ECOL):
    t.set_color(col)
ax.set_yticks([-0.4, 0, 0.4])
ax.tick_params(labelsize=5.8, pad=1.5)
ax.set_xlabel(r'delay $\tau$ (chips)', fontsize=6.4, labelpad=1.5)
ax.set_ylabel(r'wavelength $\lambda_m$ (nm)', fontsize=6.4, labelpad=0.5)
save(fig, 'panel_ridge')

# =============================== spectrum ===================================
fig = plt.figure(figsize=(1.52, 2.45))
ax = fig.add_axes([0.10, 0.125, 0.86, 0.860])
fine = np.linspace(nu[0], nu[-1], 300)
fitc = base_f + A_f * np.exp(-0.5 * ((fine - mu_f) / sig_f) ** 2)
fitc = fitc / corr[:, bins[HL]].max()
ax.plot(fitc, fine * PM / 1000.0, color=ECOL[HL], lw=1.0)
ax.plot(spec, lam, ls='none', marker='o', ms=2.6, color=ECOL[HL],
        mec='white', mew=0.35, zorder=5)
ax.axhline(mu_f * PM / 1000.0, color='0.3', ls=(0, (3, 2)), lw=0.6)
ax.text(0.02, mu_f * PM / 1000.0 + 0.035, r'$\lambda_B$', fontsize=7.2,
        color='0.15')
ax.text(0.98, YLO + 0.55 * GAP,
        r'$\Delta\lambda_B \propto \Delta T,\ \varepsilon$',
        fontsize=6.2, color='0.15', ha='right')
ax.set_xlim(-0.06, 1.10)
ax.set_ylim(YLO, YHI)
ax.set_xticks([0, 1]); ax.set_yticks([])
ax.tick_params(labelsize=5.8, pad=1.5)
ax.set_xlabel('reflectivity', fontsize=6.4, labelpad=1.5)
save(fig, 'panel_spec')

# =============================== raw-record stamp ===========================
fig = plt.figure(figsize=(1.10, 0.92))
ax = fig.add_axes([0.04, 0.05, 0.92, 0.90])
ax.imshow(rec, aspect='auto', origin='lower', cmap='Greys',
          extent=[0, NCH, -0.5, M - 0.5], interpolation='nearest',
          vmax=rec.max() * 1.05)
ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(True); sp.set_linewidth(0.5); sp.set_color('0.4')
save(fig, 'panel_stamp')

# =============================== inset: drive ===============================
OSR = 8
NSTEP, CHPS = 4, 12
stair = []
for i in range(NSTEP):
    stair.append(i * 1.0 + 0.34 * 0.5 * (1 + _MS[i * CHPS:(i + 1) * CHPS]))
stair = np.repeat(np.concatenate(stair), OSR)
ts = np.arange(len(stair)) / (OSR * CHPS)
fig = plt.figure(figsize=(1.30, 0.62))
ax = fig.add_axes([0.02, 0.06, 0.96, 0.90])
base = np.repeat(np.arange(NSTEP) * 1.0, CHPS * OSR)
ax.plot(ts, base, color=BLUE, lw=0.55, alpha=0.35)
ax.plot(ts, stair, color=BLUE, lw=0.65)
ax.set_xlim(-0.04, ts[-1] + 0.04); ax.set_ylim(-0.35, NSTEP - 0.10)
ax.axis('off')
save(fig, 'inset_drive')

# =============================== inset: echoes + sum ========================
taus = [2.2, 6.0, 10.2, 14.2]
amps = [1.00, 0.74, 0.52, 0.36]
drv2 = np.repeat(0.5 * (1 + _MS[:34]), 6)
tt2 = np.arange(len(drv2)) / 6.0
W2 = tt2 < 30
fig = plt.figure(figsize=(2.45, 1.30))
ax = fig.add_axes([0.02, 0.04, 0.96, 0.92])
offs = (5.6, 4.4, 3.2, 2.0)
ech = []
for tcau, a, col, off in zip(taus, amps, ECOL, offs):
    e = np.zeros_like(drv2); sh = int(tcau * 6)
    e[sh:] = a * drv2[:len(drv2) - sh]
    ech.append(e)
    ax.plot(tt2[W2], off + 0.85 * e[W2], color=col, lw=0.6)
sm = np.sum(ech, axis=0)
ax.plot(tt2[W2], 0.0 + 0.42 * sm[W2], color='0.1', lw=0.7)
# transmit marker and tau_2 dimension
ax.axvline(0, color='0.45', lw=0.5, ls=(0, (2, 1.5)))
y2 = offs[1] + 1.05
ax.annotate('', xy=(taus[1], y2), xytext=(0, y2),
            arrowprops=dict(arrowstyle='-|>', lw=0.6, color=ECOL[1]))
ax.text(taus[1] / 2, y2 + 0.12, r'$\tau_2$', fontsize=6.4, color=ECOL[1],
        ha='center', va='bottom')
ax.set_xlim(-1.2, 30); ax.set_ylim(-0.3, 7.1)
ax.axis('off')
save(fig, 'inset_echosum')

# =============================== inset: overlap =============================
fig = plt.figure(figsize=(1.55, 0.66))
ax = fig.add_axes([0.03, 0.07, 0.94, 0.88])
lamg = np.linspace(-40, 40, 240)
for c0, col in zip([-9.0, 3.0, -3.0, 7.0], ECOL):   # crowded around one lamB
    g = np.exp(-0.5 * ((lamg - c0) / SIG) ** 2)
    ax.plot(lamg, g, color=col, lw=0.65)
    ax.fill_between(lamg, g, color=col, alpha=0.06, lw=0)
ax.axvline(-0.5, color='0.35', lw=0.5, ls=(0, (2.5, 1.8)))
ax.set_xlim(-40, 40); ax.set_ylim(0, 1.12)
ax.axis('off')
save(fig, 'inset_overlap')

print('fitted lambda_B: %.1f pm (true %.1f pm), M=%d'
      % (mu_f * PM, nub[HL] * PM, M))
