"""s24_concept.py - data panels for the concept figure (Fig. 1), v10.

The figure is composed in TikZ (Draft_v2/fig1_concept.tex); this script
renders only the DATA panels, each as a tight standalone PDF+PNG, with STIX
(Times-like) typography matching the IEEEtran page.

Illustration uses M=16 wavelength steps so that the discrete structure is
visible everywhere: the tuning staircase is genuinely square, the raw record
has sixteen fat rows, the map stacks sixteen correlations, the spectrum panel
shows sixteen samples under a Gaussian fit. The physics model is the paper's:
m-sequence N=127, K=4 gratings at R=5%, shadowing via cumulative
transmission, leakage via the periodic autocorrelation.

Panels:
  panel_record   wavelength staircase strip + raw record (M x N*OS image)
  panel_map      the record after correlating every row
  panel_spec     highlighted column: 16 samples, Gaussian fit, lambda_B
  prof_row       1D: what the photodiode wrote during step lambda_m
  prof_corr      1D: the same row after correlation (peaks + leakage floor)
  inset_drive    staircase drive, one code period per step (square steps)
  inset_tuning   VCSEL tuning curve lambda(V), slightly bowed, step points
  inset_chirp    code chirps the laser: adiabatic level + edge transients
  inset_echoes   the K delayed, weighted copies of the code
  inset_overlap  the K overlapping reflection spectra
  inset_template one period of the transmitted sequence
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
ECOL = ['#D55E00', '#E69F00', '#009E73', '#CC79A7']

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
NCH = 127
_MS = 1.0 - 2.0 * C._mls01(7)
_MS = np.roll(_MS, -int(np.argmax(_MS == 1.0)))
ACORR = C.periodic_xcorr(_MS, _MS)

# --- model ------------------------------------------------------------------
M = 16
nu = np.linspace(-2.15 * F, 2.15 * F, M)
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
corr = (Wl @ prim).T                     # (M, NCH), full model
HL = 1
JM = np.argmin(np.abs(nu - 0.0))
spec = corr[:, bins[HL]]
A_f, mu_f, sig_f, base_f = C.gauss_fit_full(nu, spec)
mu_row = np.interp(mu_f, nu, np.arange(M))

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


# --- panel_record: staircase strip + raw record -----------------------------
fig = plt.figure(figsize=(2.28, 1.81))
axl = fig.add_axes([0.145, 0.175, 0.150, 0.805])
axb = fig.add_axes([0.325, 0.175, 0.665, 0.805])
axl.plot(lam, np.arange(M), color=BLUE, lw=0.9, drawstyle='steps-mid')
axl.plot([lam[JM]], [JM], marker='o', ms=2.4, color='0.15', zorder=5)
axl.set_ylim(-0.5, M - 0.5); axl.set_xlim(-0.55, 0.55)
axl.set_xticks([]); axl.set_yticks([0, M - 1])
axl.set_yticklabels(['1', '$M$'])
axl.tick_params(labelsize=5.5, pad=1.2)
axl.set_xlabel(r'$\lambda_m$', fontsize=6, labelpad=1)
axl.set_ylabel('wavelength step $m$', fontsize=6, labelpad=0)
axl.spines['bottom'].set_visible(False)
axb.imshow(rec, aspect='auto', origin='lower', cmap='Greys',
           extent=[0, NCH, -0.5, M - 0.5], interpolation='nearest',
           vmax=rec.max())
axb.axhline(JM, color=BLUE, lw=0.6, ls=(0, (3, 1.5)))
axb.text(4, JM + 0.75, r'step $\lambda_m$', fontsize=5.5, color=BLUE,
         bbox=dict(fc='white', ec='none', alpha=0.75, pad=0.5))
axb.set_xlim(0, NCH); axb.set_ylim(-0.5, M - 0.5)
axb.set_xticks([0, NCH]); axb.set_xticklabels(['0', '$N$'])
axb.set_yticks([])
axb.tick_params(labelsize=5.5, pad=1.2)
axb.set_xlabel('time in code period (chips)', fontsize=6, labelpad=1)
save(fig, 'panel_record')

# --- panel_map: correlated record -------------------------------------------
fig = plt.figure(figsize=(2.09, 1.81))
ax = fig.add_axes([0.06, 0.175, 0.90, 0.805])
ax.imshow(corr, aspect='auto', origin='lower', cmap='Greys',
          extent=[0, NCH, -0.5, M - 0.5], interpolation='nearest',
          vmax=corr.max() * 0.9)
for b, col in zip(bins, ECOL):
    ax.plot(b, M - 0.20, marker='v', ms=2.6, color=col, clip_on=False)
ax.axvspan(bins[HL] - 3.5, bins[HL] + 3.5, color=ECOL[HL], alpha=0.16, lw=0)
ax.axhline(JM, color=BLUE, lw=0.6, ls=(0, (3, 1.5)))
ax.text(123, JM + 0.75, 'same row', fontsize=5.5, color=BLUE, ha='right')
ax.text(bins[HL] + 5.5, 0.15, 'column $=$ grating', fontsize=5.5,
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

# --- panel_spec: highlighted column, samples + Gaussian fit -----------------
fig = plt.figure(figsize=(1.46, 2.36))
ax = fig.add_axes([0.13, 0.135, 0.80, 0.845])
sn = spec / spec.max()
fine = np.linspace(nu[0], nu[-1], 300)
fitc = (base_f + A_f * np.exp(-0.5 * ((fine - mu_f) / sig_f) ** 2))
fitc = fitc / spec.max()
ax.plot(fitc, np.interp(fine, nu, np.arange(M)), color=ECOL[HL], lw=0.9)
ax.plot(sn, np.arange(M), ls='none', marker='o', ms=2.4, color=ECOL[HL],
        mec='white', mew=0.3, zorder=5)
ax.axhline(mu_row, color='0.35', ls=(0, (3, 2)), lw=0.6)
ax.text(0.99, mu_row + 0.35, r'$\lambda_B$', fontsize=7, color='0.15',
        ha='right')
ax.annotate('Gaussian fit', xy=(fitc[205], np.interp(fine[205], nu,
            np.arange(M))), xytext=(0.98, 12.6), fontsize=5.5, color='0.15',
            ha='right',
            arrowprops=dict(arrowstyle='-', lw=0.5, color='0.6'))
ax.annotate(r'$M$ samples', xy=(sn[11], 11), xytext=(0.98, 13.9),
            fontsize=5.5, color='0.15', ha='right',
            arrowprops=dict(arrowstyle='-', lw=0.5, color='0.6'))
ax.text(0.97, 1.4, r'$\Delta\lambda_B \propto \Delta T,\ \varepsilon$',
        fontsize=6.4, color='0.15', ha='right')
ax.set_xlim(-0.05, 1.06); ax.set_ylim(-0.5, M - 0.5)
ax.set_xticks([0, 1]); ax.set_yticks([])
ax.tick_params(labelsize=5.5, pad=1.2)
ax.set_xlabel('reflectivity (norm.)', fontsize=6, labelpad=1)
save(fig, 'panel_spec')

# --- prof_row: the photodiode trace of the highlighted step -----------------
fig = plt.figure(figsize=(1.89, 0.43))
ax = fig.add_axes([0.015, 0.10, 0.97, 0.80])
ax.plot(np.arange(NCH * OS) / OS, rec[JM], color='0.25', lw=0.55)
ax.set_xlim(0, NCH); ax.set_ylim(rec[JM].min() - 0.005,
                                 rec[JM].max() + 0.005)
ax.axis('off')
save(fig, 'prof_row')

# --- prof_corr: the same row after correlation ------------------------------
fig = plt.figure(figsize=(1.89, 0.43))
ax = fig.add_axes([0.015, 0.10, 0.97, 0.80])
tr = corr[JM] / corr[JM].max()
ax.fill_between(np.arange(NCH), 0, np.abs(tr), color='0.80', lw=0)
for b, col in zip(bins, ECOL):
    ax.plot([b, b], [0, tr[b]], color=col, lw=1.2, solid_capstyle='butt')
    ax.plot([b], [tr[b]], marker='o', ms=1.9, color=col)
ax.set_xlim(0, NCH); ax.set_ylim(-0.03, 1.10)
ax.axis('off')
save(fig, 'prof_corr')

# --- inset_drive: square staircase, one code period per step ----------------
OSR = 8
NSTEP, CHPS = 4, 12
stair = []
for i in range(NSTEP):
    stair.append(i * 1.0 + 0.34 * 0.5 * (1 + _MS[i * CHPS:(i + 1) * CHPS]))
stair = np.repeat(np.concatenate(stair), OSR)
ts = np.arange(len(stair)) / (OSR * CHPS)
fig = plt.figure(figsize=(1.22, 0.60))
ax = fig.add_axes([0.02, 0.06, 0.96, 0.90])
base = np.repeat(np.arange(NSTEP) * 1.0, CHPS * OSR)
ax.plot(ts, base, color=BLUE, lw=0.55, alpha=0.35)
ax.plot(ts, stair, color=BLUE, lw=0.6)
ax.set_xlim(-0.04, ts[-1] + 0.04); ax.set_ylim(-0.35, NSTEP - 0.15)
ax.axis('off')
save(fig, 'inset_drive')

# --- inset_tuning: lambda(V), slightly bowed, with the step points ----------
fig = plt.figure(figsize=(1.22, 0.60))
ax = fig.add_axes([0.03, 0.07, 0.94, 0.88])
V = np.linspace(0, 1, 200)
lamV = V + 0.16 * V * (1 - V) * 2.0          # bowed tuning curve
ax.plot(V, lamV, color=BLUE, lw=0.7)
Vs = np.linspace(0.06, 0.94, 8)
ax.plot(Vs, Vs + 0.16 * Vs * (1 - Vs) * 2.0, ls='none', marker='o', ms=1.7,
        color='0.15')
ax.plot([0, 1], [0.0, 1.0], color='0.55', lw=0.5, ls=(0, (2, 1.5)))
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.05, 1.15)
ax.axis('off')
save(fig, 'inset_tuning')

# --- inset_chirp: adiabatic level + edge transients on one code chunk -------
fig = plt.figure(figsize=(1.22, 0.60))
ax = fig.add_axes([0.02, 0.06, 0.96, 0.90])
tt = np.linspace(0, 4, 800)
P = ((tt % 2) < 1).astype(float)              # chips: on-off-on-off
dnu = 0.55 * P.copy()
edges = np.where(np.abs(np.diff(P)) > 0)[0]
for e in edges:
    sgn = np.sign(P[e + 1] - P[e])
    tail = np.exp(-(tt[e:] - tt[e]) / 0.09)
    dnu[e:] += 0.45 * sgn * tail
ax.fill_between(tt, -0.35, 1.25, where=P > 0.5, color='0.92', lw=0)
ax.plot(tt, dnu, color=BLUE, lw=0.6)
ax.axhline(0.55, color='0.55', lw=0.45, ls=(0, (2, 1.5)))
ax.set_xlim(0, 4); ax.set_ylim(-0.55, 1.25)
ax.axis('off')
save(fig, 'inset_chirp')

# --- inset_echoes -----------------------------------------------------------
taus = [3.0, 7.5, 12.0, 16.5]
drv2 = np.repeat(0.5 * (1 + _MS[:40]), 6)
tt2 = np.arange(len(drv2)) / 6.0
W2 = tt2 < 34
fig = plt.figure(figsize=(1.22, 0.60))
ax = fig.add_axes([0.02, 0.06, 0.96, 0.90])
for tcau, a, col, off in zip(taus, [1.0, 0.80, 0.62, 0.48], ECOL,
                             (3.6, 2.4, 1.2, 0.0)):
    e = np.zeros_like(drv2); sh = int(tcau * 6)
    e[sh:] = a * drv2[:len(drv2) - sh]
    ax.plot(tt2[W2], off + e[W2], color=col, lw=0.5)
ax.set_xlim(0, 34); ax.set_ylim(-0.3, 5.0)
ax.axis('off')
save(fig, 'inset_echoes')

# --- inset_overlap ----------------------------------------------------------
fig = plt.figure(figsize=(1.22, 0.60))
ax = fig.add_axes([0.03, 0.07, 0.94, 0.88])
lamg = np.linspace(-60, 60, 200)
for c0, col in zip(nub, ECOL):
    g = np.exp(-0.5 * ((lamg - c0) / SIG) ** 2)
    ax.plot(lamg, g, color=col, lw=0.55)
    ax.fill_between(lamg, g, color=col, alpha=0.10, lw=0)
ax.set_xlim(-60, 60); ax.set_ylim(0, 1.08)
ax.axis('off')
save(fig, 'inset_overlap')

# --- inset_template ---------------------------------------------------------
fig = plt.figure(figsize=(0.85, 0.34))
ax = fig.add_axes([0.02, 0.10, 0.96, 0.80])
ax.plot(np.arange(len(code01)) / OS, code01, color='0.3', lw=0.55)
ax.set_xlim(0, 22); ax.set_ylim(-0.25, 1.25)
ax.axis('off')
save(fig, 'inset_template')

print('fitted lambda_B: %.1f pm (true %.1f pm), M=%d steps'
      % (mu_f * PM, nub[HL] * PM, M))
