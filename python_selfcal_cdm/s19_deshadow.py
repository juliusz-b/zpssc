"""s19_deshadow.py - the dominant array error, and the recursion that removes it.

Spectral shadowing is the term that binds when the gratings are strong, and it is
also the one term that can be undone rather than merely tolerated, because it is
deterministic once the gratings in front are known. Three panels.

  (a) The recursion. The first grating is illuminated by the full source, so its
      line is fitted directly. From the fitted amplitude, centre and width the
      transmission it imposes is reconstructed and divided out of every grating
      behind it. The second grating is then clean, and the procedure repeats.
      Nothing is needed beyond the order of the gratings, which the delay bins
      already give.

  (b) One grating, three curves: the true line, what the array actually returns
      once three strong gratings sit in front of it, and what the recursion
      recovers. The fitted peak moves back to where it belongs.

  (c) How far the recursion gets. It works while each line can still be fitted
      cleanly, which holds to roughly twenty gratings; past that, each stage
      feeds its own error into the next and the gain disappears. Shadowing, code
      leakage and detector noise are included, so the recursion is judged against
      a floor it cannot remove; ghosts are left out because they have their own
      figure.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import warnings; warnings.filterwarnings('ignore')
import common as C

plt.rcParams.update({'font.size': 8.5})

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
M = 96
DETUNE = 25.0
nu = np.linspace(-2.6 * F, 2.6 * F, M)
SIGMA_N = 1.1e-6


def shadowed_array(K, R, rng):
    """Per-grating readout with cumulative upstream transmission."""
    nub = rng.uniform(-DETUNE, DETUNE, size=K)
    shapes = np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
    tcum = np.ones((K, M))
    for k in range(1, K):
        tcum[k] = tcum[k - 1] * (1.0 - R * shapes[k - 1]) ** 2
    return nub, R * shapes, R * shapes * tcum


def peel(S, floor=0.05):
    """Sequential deshadowing: fit, reconstruct the transmission, divide it out."""
    K = S.shape[0]
    tcorr = np.ones(M)
    out = np.empty(K)
    corrected = np.empty_like(S)
    for k in range(K):
        Sk = S[k] / np.maximum(tcorr, floor)
        corrected[k] = Sk
        a, mu, sg, _ = C.gauss_fit_full(nu, Sk)
        out[k] = mu
        line = np.clip(a, 0.0, 0.99) * np.exp(-0.5 * ((nu - mu) / max(sg, 1e-3)) ** 2)
        tcorr = tcorr * (1.0 - line) ** 2
    return out, corrected


NCH = 127
_MS = 1.0 - 2.0 * C._mls01(7)
ACORR = C.periodic_xcorr(_MS, _MS)


def rms_error(K, R, seed, corrected=False, ntrials=6):
    """Shadowing plus code leakage and noise, so the recursion is judged against
    a realistic floor rather than against an empty one. Ghosts are left out here;
    they have their own figure."""
    errs = []
    for t in range(ntrials):
        rng = np.random.default_rng(seed + t)
        nub, clean, shad = shadowed_array(K, R, rng)
        bins = np.sort(rng.choice(np.arange(1, NCH), size=K, replace=False))
        W = ACORR[(bins[:, None] - bins[None, :]) % NCH]
        np.fill_diagonal(W, 0.0)
        S = shad + W @ shad + rng.normal(0, SIGMA_N, shad.shape)
        if corrected:
            mu, _ = peel(S)
        else:
            mu = np.array([C.gauss_fit_peak(nu, S[k]) for k in range(K)])
        errs.append(np.sqrt(np.mean(((mu - nub) * PM) ** 2)))
    return float(np.mean(errs))


# --- data for panel (b): the fourth grating behind three strong ones ---------
rng_b = np.random.default_rng(12)
K_B, R_B = 4, 0.20
nub_b, clean_b, shad_b = shadowed_array(K_B, R_B, rng_b)
_, corrected_b = peel(shad_b)
kk = K_B - 1
p_true = C.gauss_fit_peak(nu, clean_b[kk]) * PM
p_shad = C.gauss_fit_peak(nu, shad_b[kk]) * PM
p_corr = C.gauss_fit_peak(nu, corrected_b[kk]) * PM

# --- data for panel (c) ------------------------------------------------------
Ks = np.array([3, 4, 8, 16, 24, 32, 48])
raw = np.array([rms_error(K, 0.10, 200 + K) for K in Ks])
fixed = np.array([rms_error(K, 0.10, 200 + K, corrected=True) for K in Ks])

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(10.6, 3.5))

# --- (a) the recursion, drawn ----------------------------------------------
axa = ax[0]
axa.axis('off'); axa.set_xlim(0, 10); axa.set_ylim(0, 10.6)
axa.plot([0.6, 9.4], [8.9, 8.9], color='#555', lw=1.5)
gz = [1.8, 4.0, 6.2, 8.4]
for i, z in enumerate(gz):
    axa.add_patch(Rectangle((z - 0.22, 8.62), 0.44, 0.56, fc='#c0392b',
                            ec='#7b241c', hatch='///', lw=0.8))
    axa.text(z, 9.45, '%d' % (i + 1), ha='center', fontsize=8)
axa.annotate('', xy=(9.4, 9.9), xytext=(0.6, 9.9),
             arrowprops=dict(arrowstyle='-|>', lw=0.9, color='#555'))
axa.text(0.6, 9.95, 'light', ha='left', fontsize=6.8, color='#555')


def step(y, text, fc='#eaf1fb'):
    axa.add_patch(FancyBboxPatch((0.5, y), 9.0, 1.05,
                                 boxstyle='round,pad=0.05,rounding_size=0.12',
                                 fc=fc, ec='#1f4e79', lw=1.0))
    axa.text(5.0, y + 0.52, text, ha='center', va='center', fontsize=7.2)


step(6.65, 'fit grating 1: it sees the full source', )
step(4.85, r'reconstruct its transmission $[1-R_1(\lambda)]^2$', fc='#fdf3e7')
step(3.05, 'divide it out of gratings 2 ... K')
step(1.25, 'repeat with grating 2, then 3, ...', fc='#eafaf0')
for y in (6.65, 4.85, 3.05):
    axa.add_patch(FancyArrowPatch((5.0, y), (5.0, y - 0.72), arrowstyle='-|>',
                                  mutation_scale=10, lw=1.0, color='#333'))
axa.add_patch(FancyArrowPatch((5.0, 8.55), (5.0, 7.75), arrowstyle='-|>',
                              mutation_scale=10, lw=1.0, color='#333'))
axa.text(0.6, 0.45, 'needs only the order of the gratings, which the delay '
                    'bins already give', fontsize=6.6, color='0.35')
axa.set_title('(a) Sequential deshadowing', fontsize=9)

# --- (b) one grating, before and after --------------------------------------
axb = ax[1]
lam = nu * PM / 1000.0
axb.plot(lam, clean_b[kk] / clean_b[kk].max(), color='#2980b9', lw=2.2,
         label='true line')
axb.plot(lam, shad_b[kk] / shad_b[kk].max(), color='#c0392b', lw=1.4,
         label='as measured (shadowed)')
axb.plot(lam, corrected_b[kk] / corrected_b[kk].max(), color='#f39c12', lw=1.4,
         ls='--', label='after deshadowing')
axb.axvline(p_true / 1000.0, color='#2980b9', ls=':', lw=0.8)
axb.axvline(p_shad / 1000.0, color='#c0392b', ls=':', lw=0.8)
axb.set_xlim(-0.42, 0.14); axb.set_ylim(0, 1.32)
axb.set_xlabel('wavelength offset [nm]'); axb.set_ylabel('normalised readout')
axb.set_title('(b) Fourth grating behind three at R = 20%', fontsize=9)
axb.legend(fontsize=6.6, loc='upper left')
axb.text(0.97, 0.10, 'peak error %.0f pm  ->  %.1f pm'
         % (abs(p_shad - p_true), abs(p_corr - p_true)),
         transform=axb.transAxes, ha='right', fontsize=7.2)

# --- (c) how far it gets -----------------------------------------------------
axc = ax[2]
axc.semilogy(Ks, raw, 'o-', color='#c0392b', label='uncorrected')
axc.semilogy(Ks, fixed, 's-', color='#f39c12', label='after deshadowing')
axc.axhline(10.0, color='0.3', ls='--', lw=0.8)
axc.text(0.97, 0.06, '10 pm target', transform=axc.transAxes, fontsize=7,
         color='0.3', ha='right')
axc.set_xlabel('gratings on the fiber, K')
axc.set_ylabel('RMS Bragg error [pm]')
axc.set_title('(c) Gain and its limit, R = 10%', fontsize=9)
axc.legend(fontsize=7, loc='upper left'); axc.grid(True, which='both', alpha=0.25)

fig.tight_layout()
fig.savefig('figs/fig_s19_deshadow.png', dpi=150, bbox_inches='tight')

print('panel (b): K=%d, R=%.2f, grating %d' % (K_B, R_B, kk + 1))
print('  true %.1f pm, shadowed %.1f pm, deshadowed %.1f pm'
      % (p_true, p_shad, p_corr))
print('--- (c) shadowing with code leakage and noise, R = 10% ---')
print('     K   uncorrected [pm]   deshadowed [pm]   gain')
for K, a, b in zip(Ks, raw, fixed):
    print('  %4d   %16.2f   %15.2f   %5.1fx' % (K, a, b, a / max(b, 1e-9)))
print('saved figs/fig_s19_deshadow.png')
