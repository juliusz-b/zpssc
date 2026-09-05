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
import figstyle as FS
FS.apply()


PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
M = 96
DETUNE = 25.0
nu = np.linspace(-2.6 * F, 2.6 * F, M)
SIGMA_N = 1.1e-6


def shadowed_array(K, R, rng, nub=None):
    """Per-grating readout with cumulative upstream transmission."""
    if nub is None:
        nub = rng.uniform(-DETUNE, DETUNE, size=K)
    nub = np.asarray(nub, float)
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
# panel (a): odstrojenia na sztywno (GHz), tak zeby kazda z trzech siatek z przodu
# widocznie zacieniala czwarta i kazdy krok korekcji byl widoczny
nub_b, clean_b, shad_b = shadowed_array(K_B, R_B, rng_b, nub=[-7.0, 8.0, 4.0, 0.0])
_, corrected_b = peel(shad_b)
kk = K_B - 1


def peel_steps(S, kk, floor=0.05):
    """Czesciowe korekcje siatki kk: S_kk / (T1), / (T1 T2), ... jak w peel()."""
    tcorr = np.ones(M)
    steps = []
    for k in range(kk):
        Sk = S[k] / np.maximum(tcorr, floor)
        a, mu, sg, _ = C.gauss_fit_full(nu, Sk)
        line = np.clip(a, 0.0, 0.99) * np.exp(-0.5 * ((nu - mu) / max(sg, 1e-3)) ** 2)
        tcorr = tcorr * (1.0 - line) ** 2
        steps.append(S[kk] / np.maximum(tcorr, floor))
    return steps


steps_b = peel_steps(shad_b, kk)
p_steps = [C.gauss_fit_peak(nu, st) * PM for st in steps_b]
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
# two panels: the recursion itself lives in eqs. (19)-(20) of the paper
fig, ax = plt.subplots(2, 1, figsize=(2.7, 3.6))

# --- (a) one grating, the correction step by step ---------------------------
axb = ax[0]
lam = nu * PM / 1000.0
nrm = clean_b[kk].max()
axb.plot(lam, clean_b[kk] / nrm, color='0.55', lw=2.6, alpha=0.6, label='$R_4$')
axb.plot(lam, shad_b[kk] / nrm, color='#D55E00', lw=1.4, label='$S_4$')
step_cols = ['#E69F00', '#56B4E9']
step_lab = ['$S_4/\\widehat T_1$', '$S_4/\\widehat T_1\\widehat T_2$']
for st, col, lab in zip(steps_b[:-1], step_cols, step_lab):
    axb.plot(lam, st / nrm, color=col, lw=1.0, label=lab)
axb.plot(lam, corrected_b[kk] / nrm, color='#009E73', lw=1.4, ls='--', label='$\\widehat S_4$')
axb.axvline(p_true / 1000.0, color='0.3', ls=':', lw=0.8)
axb.axvline(p_shad / 1000.0, color='#D55E00', ls=':', lw=0.8)
FS.dim_gap(axb, p_true / 1000.0, p_shad / 1000.0, 1.13,
           '%.0f pm' % abs(p_shad - p_true), color='#D55E00',
           tail=0.055, side='right')
axb.set_xlim(-0.45, 0.45); axb.set_ylim(0, 1.62)
axb.set_xlabel('wavelength offset [nm]'); axb.set_ylabel('readout / peak of $R_4$')
FS.letter(axb, 'a')
axb.legend(fontsize=5.5, loc='upper left', ncol=1, frameon=True,
           handlelength=1.4, columnspacing=0.7, labelspacing=0.2)

# --- (b) how far it gets -----------------------------------------------------
axc = ax[1]
axc.semilogy(Ks, raw, 'o-', color='#D55E00', label='uncorrected $S_k$')
axc.semilogy(Ks, fixed, 's-', color='#009E73', label='corrected $\\widehat S_k$')
axc.axhline(10.0, color='0.3', ls='--', lw=0.8, label='10 pm target')
axc.set_xlabel('gratings on the fiber, $K$')
axc.set_ylabel('RMS $\\delta\\lambda_k$ [pm]')
FS.letter(axc, 'b')
axc.legend(fontsize=5.5, loc='lower right',
           ncol=1, frameon=True, handlelength=1.4, columnspacing=0.65)
axc.grid(False, which='both', alpha=0.25)

fig.subplots_adjust(left=0.19, right=0.98, top=0.95, bottom=0.10, hspace=0.42)
fig.savefig('figs/fig_s19_deshadow.png', dpi=150, bbox_inches='tight')
fig.savefig('figs/fig_s19_deshadow.pdf', bbox_inches='tight')

print('panel (b): K=%d, R=%.2f, grating %d' % (K_B, R_B, kk + 1))
print('  true %.1f pm, shadowed %.1f pm, deshadowed %.1f pm'
      % (p_true, p_shad, p_corr))
print('--- (c) shadowing with code leakage and noise, R = 10% ---')
print('grating 4 behind three at R=20%%: true %.1f, shadowed %.1f, after T1 %.1f, after T1T2 %.1f, corrected %.1f pm'
      % (p_true, p_shad, p_steps[0], p_steps[1], p_corr))
print('     K   uncorrected [pm]   deshadowed [pm]   gain')
for K, a, b in zip(Ks, raw, fixed):
    print('  %4d   %16.2f   %15.2f   %5.1fx' % (K, a, b, a / max(b, 1e-9)))
print('saved figs/fig_s19_deshadow.png')
