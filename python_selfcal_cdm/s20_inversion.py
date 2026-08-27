"""s20_inversion.py - how much does the correction need to know, and what else is
in the despread map?

Prompted by the obvious objection to the deshadowing of s19: it rebuilds each
grating's transmission from a Gaussian fit, so it looks as if it needs the
lineshape. It does not, and two further questions follow from asking that one.

  1. SHAPE DEPENDENCE. The transmission a grating imposes is a function of its
     reflection spectrum, and the despread readout already IS that spectrum up to
     a scale factor. The correction can therefore be built point by point from
     the measurement,

         T_k(lambda) = [1 - g * S_k(lambda)]^2 ,

     with g the only thing that has to be known: the gain that converts readout
     units into reflectance. The parametric variant of s19 is a denoised special
     case, and it is biased when the assumed shape is not the true one. Both are
     compared on Gaussian and on tanh^2 gratings, and the non-parametric one is
     then stressed with a wrong g.

  2. GHOSTS ARE IDENTIFIABLE, NOT JUST HARMFUL. A third-order ghost is the
     product of three grating spectra, so its line is narrower by sqrt(3) and its
     amplitude scales as R^3. Either feature separates it from a real grating in
     the model; the delay-triple relation is a third, independent check. This is
     the natural place to ask whether a learned classifier is needed, and the
     answer here is that it is not.

  3. TWO GRATINGS IN ONE DELAY BIN. The chip rate is the expensive part of the
     receiver, and it is set by the requirement that no two gratings share a
     delay bin. That requirement can be bought off with wavelength instead: when
     two co-binned gratings differ in Bragg wavelength by at least half a
     linewidth, a two-component fit separates them. Below that the fit breaks,
     and that is exactly the regime the learned demodulators in the FBG
     literature were built for.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
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
NCH = 127
_MS = 1.0 - 2.0 * C._mls01(7)
ACORR = C.periodic_xcorr(_MS, _MS)


def build(K, R, rng, shape='gauss'):
    """True detunings, reflection spectra and the shadowed, leaky readout."""
    nub = rng.uniform(-DETUNE, DETUNE, size=K)
    if shape == 'gauss':
        spec = np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
    else:
        spec = np.array([C.fbg_tanh(nu, x, F, n_side=0.0) for x in nub])
    refl = R * spec
    tcum = np.ones((K, M))
    for k in range(1, K):
        tcum[k] = tcum[k - 1] * (1.0 - refl[k - 1]) ** 2
    bins = np.sort(rng.choice(np.arange(1, NCH), size=K, replace=False))
    S = refl * tcum
    W = ACORR[(bins[:, None] - bins[None, :]) % NCH]
    np.fill_diagonal(W, 0.0)
    S = S + W @ S + rng.normal(0, SIGMA_N, S.shape)
    return nub, refl, S, bins


def peel_parametric(S, floor=0.05):
    """s19 variant: rebuild each line from a Gaussian fit."""
    K = S.shape[0]
    t = np.ones(M)
    out = np.empty(K)
    for k in range(K):
        Sk = S[k] / np.maximum(t, floor)
        a, mu, sg, _ = C.gauss_fit_full(nu, Sk)
        out[k] = mu
        line = np.clip(a, 0.0, 0.99) * np.exp(-0.5 * ((nu - mu) / max(sg, 1e-3)) ** 2)
        t = t * (1.0 - line) ** 2
    return out


def peel_nonparametric(S, gain=1.0, floor=0.05):
    """Shape-free variant: the measured spectrum supplies the transmission."""
    K = S.shape[0]
    t = np.ones(M)
    out = np.empty(K)
    for k in range(K):
        Sk = S[k] / np.maximum(t, floor)
        out[k] = C.gauss_fit_peak(nu, Sk)
        t = t * (1.0 - np.clip(gain * Sk, 0.0, 0.99)) ** 2
    return out


def rms(K, R, seed, method, shape='gauss', gain=1.0, ntrials=6):
    errs = []
    for t in range(ntrials):
        rng = np.random.default_rng(seed + t)
        nub, refl, S, bins = build(K, R, rng, shape=shape)
        if method == 'none':
            mu = np.array([C.gauss_fit_peak(nu, S[k]) for k in range(K)])
        elif method == 'param':
            mu = peel_parametric(S)
        else:
            mu = peel_nonparametric(S, gain=gain)
        errs.append(np.sqrt(np.mean(((mu - nub) * PM) ** 2)))
    return float(np.mean(errs))


# ---------------------------------------------------------------------------
# 1. lineshape dependence and gain tolerance
# ---------------------------------------------------------------------------
K1, R1 = 8, 0.10
shape_tab = []
for shape in ('gauss', 'tanh'):
    shape_tab.append((shape, [rms(K1, R1, 400, 'none', shape=shape),
                              rms(K1, R1, 400, 'param', shape=shape),
                              rms(K1, R1, 400, 'nonparam', shape=shape)]))
gains = np.array([0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.6, 2.0])
gain_curve = np.array([rms(K1, R1, 400, 'nonparam', shape='tanh', gain=g)
                       for g in gains])

# ---------------------------------------------------------------------------
# 2. ghost against grating, from the map alone
# ---------------------------------------------------------------------------
def ghost_features(ntrials=40, K=10, R=0.20, seed=5, noise=1.0e-3):
    """Fitted width and amplitude of real gratings and of third-order ghosts,
    with a manufacturing spread on both linewidth and reflectance so that the
    two classes are compared as they would be measured, not as idealised."""
    wr, ar, wg, ag = [], [], [], []
    for t in range(ntrials):
        rng = np.random.default_rng(seed + t)
        nub = rng.uniform(-DETUNE, DETUNE, size=K)
        sig_k = SIG * rng.uniform(0.90, 1.10, size=K)     # +-10% linewidth spread
        R_k = R * rng.uniform(0.80, 1.20, size=K)         # +-20% reflectance spread
        for k in range(K):
            line = R_k[k] * np.exp(-0.5 * ((nu - nub[k]) / sig_k[k]) ** 2)
            a, mu, sg, _ = C.gauss_fit_full(nu, line + rng.normal(0, noise, M))
            wr.append(sg / SIG); ar.append(a / R)
        for _ in range(K):
            i, j, l = rng.choice(K, 3, replace=True)
            s2 = 1.0 / (1.0 / sig_k[i] ** 2 + 1.0 / sig_k[j] ** 2 + 1.0 / sig_k[l] ** 2)
            cbar = s2 * (nub[i] / sig_k[i] ** 2 + nub[j] / sig_k[j] ** 2 +
                         nub[l] / sig_k[l] ** 2)
            expo = (((nub[i] - cbar) / sig_k[i]) ** 2 + ((nub[j] - cbar) / sig_k[j]) ** 2
                    + ((nub[l] - cbar) / sig_k[l]) ** 2)
            amp = R_k[i] * R_k[j] * R_k[l] * np.exp(-0.5 * expo)
            line = amp * np.exp(-0.5 * ((nu - cbar) / np.sqrt(s2)) ** 2)
            a, mu, sg, _ = C.gauss_fit_full(nu, line + rng.normal(0, noise, M))
            wg.append(sg / SIG); ag.append(a / R)
    return np.array(wr), np.array(ar), np.array(wg), np.array(ag)


NOISE = 1.0e-3
R_GH = 0.20
wr, ar, wg, ag = ghost_features(noise=NOISE, R=R_GH)
# only peaks that are detected at all can be classified: gate at five times the
# noise, below which a ghost is not a peak but part of the floor
gate = 5 * NOISE / R_GH
det_r, det_g = ar > gate, ag > gate
wr, ar = wr[det_r], ar[det_r]
wg, ag = wg[det_g], ag[det_g]
frac_det = float(np.mean(det_g))
thr_w = 0.5 * (np.median(wr) + np.median(wg))
acc_w = 0.5 * (np.mean(wr > thr_w) + np.mean(wg < thr_w))

# ---------------------------------------------------------------------------
# 3. two gratings sharing one delay bin
# ---------------------------------------------------------------------------
nu_wide = np.linspace(-3.2 * F, 3.2 * F, M)


def _two(x, a1, m1, s1, a2, m2, s2, c):
    return (a1 * np.exp(-0.5 * ((x - m1) / s1) ** 2) +
            a2 * np.exp(-0.5 * ((x - m2) / s2) ** 2) + c)


def co_binned(sep_frac, ntrials=120, seed=7):
    """RMS error and success rate of a two-component fit on one shared bin."""
    rng = np.random.default_rng(seed)
    sep = sep_frac * F
    e1, e2, ok = [], [], 0
    for _ in range(ntrials):
        c0 = rng.uniform(-15, 15)
        m1, m2 = c0 - sep / 2, c0 + sep / 2
        a1, a2 = 0.10, 0.10 * rng.uniform(0.7, 1.3)
        y = (a1 * np.exp(-0.5 * ((nu_wide - m1) / SIG) ** 2) +
             a2 * np.exp(-0.5 * ((nu_wide - m2) / SIG) ** 2) +
             rng.normal(0, 3e-4, M))
        e1.append((C.gauss_fit_peak(nu_wide, y) - m1) * PM)
        try:
            p, _ = curve_fit(_two, nu_wide, y,
                             p0=[a1, m1 - 2, SIG, a2, m2 + 2, SIG, 0.0], maxfev=40000)
            mm = sorted([p[1], p[4]])
            if abs(mm[0] - m1) < 0.6 * sep and abs(mm[1] - m2) < 0.6 * sep:
                ok += 1
                e2.append((mm[0] - m1) * PM)
        except Exception:
            pass
    e1 = np.array(e1)
    e2 = np.array(e2) if len(e2) else np.array([np.nan])
    return (float(np.sqrt(np.mean(e1 ** 2))), float(np.sqrt(np.nanmean(e2 ** 2))),
            ok / ntrials)


seps = np.array([0.15, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0])
co = np.array([co_binned(s) for s in seps])

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(7.1, 2.35))

labels = ['uncorrected', 'parametric\n(Gaussian fit)', 'shape-free\n(measured line)']
x = np.arange(3)
w = 0.36
ax[0].bar(x - w / 2, shape_tab[0][1], w, color='#2980b9', label='Gaussian gratings')
ax[0].bar(x + w / 2, shape_tab[1][1], w, color='#c0392b', label=r'tanh$^2$ gratings')
for xi, v in zip(x - w / 2, shape_tab[0][1]):
    ax[0].text(xi, v * 1.12, '%.1f' % v, ha='center', va='bottom',
               fontsize=6.6, color='#2980b9')
for xi, v in zip(x + w / 2, shape_tab[1][1]):
    ax[0].text(xi, v * 1.12, '%.1f' % v, ha='center', va='bottom',
               fontsize=6.6, color='#c0392b')
ax[0].set_xticks(x); ax[0].set_xticklabels(labels, fontsize=7)
ax[0].set_yscale('log'); ax[0].set_ylim(0.5, 60)
ax[0].set_ylabel('RMS Bragg error [pm]')
ax[0].set_title('(a) Does deshadowing need the lineshape?', fontsize=9)
ax[0].legend(fontsize=7, loc='upper right')
ax[0].grid(True, axis='y', which='both', alpha=0.25)

ax[1].scatter(wr, ar, s=8, color='#2980b9', alpha=0.5, label='real gratings')
ax[1].scatter(wg, ag, s=8, color='#c0392b', alpha=0.5, label='third-order ghosts')
ax[1].axvline(1 / np.sqrt(3), color='0.4', ls=':', lw=0.9)
ax[1].text(1 / np.sqrt(3) + 0.015, 0.45, r'$1/\sqrt{3}$', fontsize=7, color='0.4')
ax[1].set_yscale('log')
ax[1].set_xlim(0.35, 1.25); ax[1].set_ylim(0.015, 4.0)
ax[1].set_xlabel('fitted width / grating width')
ax[1].set_ylabel('fitted amplitude / R')
ax[1].set_title('(b) A ghost looks different', fontsize=9)
ax[1].legend(fontsize=7, loc='lower right'); ax[1].grid(True, which='both', alpha=0.25)

ax[2].semilogy(seps, co[:, 0], 'o-', color='#c0392b', label='single-peak fit')
ax[2].semilogy(seps, co[:, 1], 's-', color='#28b463', label='two-component fit')
ax[2].axhline(10.0, color='0.3', ls='--', lw=0.8)
ax[2].axvspan(seps[0], 0.5, color='0.85', alpha=0.6)
ax[2].text(0.33, 40, 'fit fails\nhere', fontsize=5.8, color='0.35', ha='center')
ax[2].set_xlabel('wavelength separation / linewidth')
ax[2].set_ylabel('RMS Bragg error [pm]')
ax[2].set_title('(c) Two gratings in one delay bin', fontsize=9)
ax[2].legend(fontsize=7, loc='lower left'); ax[2].grid(True, which='both', alpha=0.25)

fig.tight_layout()
fig.savefig('figs/fig_s20_inversion.png', dpi=150, bbox_inches='tight')

# ---------------------------------------------------------------------------
# printed results
# ---------------------------------------------------------------------------
print('=== 1. lineshape dependence, K=%d, R=%.0f%% ===' % (K1, R1 * 100))
print('   truth        uncorrected   parametric   shape-free')
for shape, row in shape_tab:
    print('   %-10s   %11.2f   %10.2f   %10.2f' % (shape, row[0], row[1], row[2]))
print('   gain tolerance (tanh^2 truth, shape-free):')
for g, v in zip(gains, gain_curve):
    print('     g = %.2f -> %6.2f pm%s' % (g, v, '   <- true' if g == 1.0 else ''))
ok = gains[gain_curve < shape_tab[1][1][0]]
print('   still an improvement for g between %.1f and %.1f' % (ok.min(), ok.max()))
print()
print('=== 2. ghost against grating ===')
print('   only peaks above five times the noise are classified;'
      ' %.0f%% of ghosts reach that' % (100 * frac_det))
print('   fitted width / true width: gratings %.3f +- %.3f, ghosts %.3f +- %.3f'
      % (wr.mean(), wr.std(), wg.mean(), wg.std()))
print('   predicted ghost width ratio 1/sqrt(3) = %.3f' % (1 / np.sqrt(3)))
print('   width threshold %.3f -> balanced accuracy %.3f' % (thr_w, acc_w))
print()
print('=== 3. two gratings in one delay bin ===')
print('   sep/FWHM   sep[pm]   single-peak [pm]   two-component [pm]   success')
for s, row in zip(seps, co):
    print('   %8.2f   %7.0f   %16.1f   %18.2f   %6.0f%%'
          % (s, s * F * PM, row[0], row[1], 100 * row[2]))
print('   above half a linewidth the pair is resolved; below it the parametric')
print('   fit breaks, which is where a learned demodulator would earn its place')
print('saved figs/fig_s20_inversion.png')
