"""Is the sum of the absolute single-mechanism errors an upper bound on the full-model error?

Same model as s12_capacity (shadowing, third-order ghosts, m-sequence leakage, local Gaussian
fit, no noise), evaluated on identical random layouts in four configurations: full, shadowing
only, ghosts only, leakage only. For every grating the triangle sum
    T_k = |e_shadow,k| + |e_ghost,k| + |e_leak,k|
is compared with the full-model error |e_full,k|, and with the analytic worst case
    B_k = 1.5 sum_j |Rule A(Delta_jk)| + 0.86 a_g,k sigma + 0.86 (K-1) sigma / (N T_k),
(the factor 1.5 covers the local Gaussian fit, which responds more strongly than the full-spectrum
fit behind Rule A: without it up to 23 % of the gratings exceed the bound at R = 10 %, K = 4)
with T_k the cumulative two-pass transmission at the centre of grating k (leakage comes from the
strong front gratings, so it must be compared with the attenuated return of grating k).
To first order the shifts add, so T_k and B_k should bound |e_full,k| as long as every
mechanism is a small perturbation (R below R_c, no peak splitting).

Output: figs/fig_s55_error_bound.{pdf,png}, numbers printed for the paper.
"""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS
FS.apply()

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
NBITS = 7
N_CHIPS = 2 ** NBITS - 1
M_SLOTS = 64
BAND_GHZ = 2.6 * F
DETUNE_GHZ = 25.0
nu = np.linspace(-BAND_GHZ, BAND_GHZ, M_SLOTS)
_MSEQ = 1.0 - 2.0 * C._mls01(NBITS)
ACORR = C.periodic_xcorr(_MSEQ, _MSEQ)
SIG_PM = SIG * PM
CMAX = np.sqrt(2.0) * np.exp(-0.5)          # max of |u exp(-u^2/4)| over u = Delta/sigma, 0.858


def layout(K, rng):
    return np.sort(rng.choice(np.arange(1, N_CHIPS), size=K, replace=False))


def spectra(K, R, rng):
    nub = rng.uniform(-DETUNE_GHZ, DETUNE_GHZ, size=K)
    shapes = np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
    tcum = np.ones((K, M_SLOTS))
    trans = (1.0 - R * shapes) ** 2
    for k in range(1, K):
        tcum[k] = tcum[k - 1] * trans[k - 1]
    return nub, shapes, tcum


def ghosts(K, nub, R, bins, tcum):
    """Ghost spectra per target bin [K, M] and the ghost-to-direct amplitude ratio per grating."""
    idx = np.arange(K)
    a, b, c = np.meshgrid(idx, idx, idx, indexing='ij')
    m = (b < a) & (b < c)
    a, b, c = a[m], b[m], c[m]
    gbin = bins[a] - bins[b] + bins[c]
    pos = np.clip(np.searchsorted(bins, gbin), 0, K - 1)
    hit = bins[pos] == gbin
    out = np.zeros((K, M_SLOTS))
    ratio = np.zeros(K)
    if not hit.any():
        return out, ratio
    a, b, c, target = a[hit], b[hit], c[hit], pos[hit]
    cbar = (nub[a] + nub[b] + nub[c]) / 3.0
    spread = (nub[a] - cbar) ** 2 + (nub[b] - cbar) ** 2 + (nub[c] - cbar) ** 2
    amp = R ** 3 * np.exp(-0.5 * spread / SIG ** 2)
    up = np.minimum(np.minimum(a, b), c)
    sig_g = SIG / np.sqrt(3.0)
    for i in range(len(amp)):
        g = amp[i] * tcum[up[i]] * np.exp(-0.5 * ((nu - cbar[i]) / sig_g) ** 2)
        out[target[i]] += g
        k = target[i]
        direct_peak = R * np.interp(nub[k], nu, tcum[k])
        ratio[k] += g.max() / direct_peak
    return out, ratio


def errors(S, nub):
    return np.array([(C.gauss_fit_peak(nu, S[k]) - nub[k]) * PM for k in range(len(nub))])


def rule_a_sum(nub, R):
    """Sum over upstream gratings of |first-order shadowing shift| [pm], equal widths."""
    K = len(nub)
    out = np.zeros(K)
    for k in range(K):
        d = (nub[:k] - nub[k]) / SIG
        out[k] = np.sum(np.abs((4.0 / 3.0) * np.sqrt(2.0 / 3.0) * R * d * np.exp(-d * d / 3.0))) * SIG_PM
    return out


def one_layout(K, R, rng):
    bins = layout(K, rng)
    nub, shapes, tcum = spectra(K, R, rng)
    ones = np.ones_like(tcum)
    W = ACORR[(bins[:, None] - bins[None, :]) % N_CHIPS]
    np.fill_diagonal(W, 0.0)
    direct_sh = R * shapes * tcum
    direct_0 = R * shapes
    G_sh, ratio = ghosts(K, nub, R, bins, tcum)
    G_0, _ = ghosts(K, nub, R, bins, ones)
    e_full = errors(direct_sh + W @ direct_sh + G_sh, nub)
    e_sh = errors(direct_sh, nub)
    e_gh = errors(direct_0 + G_0, nub)
    e_lk = errors(direct_0 + W @ direct_0, nub)
    tri = np.abs(e_sh) + np.abs(e_gh) + np.abs(e_lk)
    Tk = np.array([np.interp(nub[k], nu, tcum[k]) for k in range(K)])
    ana = 1.5 * rule_a_sum(nub, R) + CMAX * ratio * SIG_PM + CMAX * (K - 1) * SIG_PM / (N_CHIPS * Tk)
    return np.abs(e_full), tri, ana


if __name__ == '__main__':
    Ks = [4, 8, 16, 24, 32, 48, 64]
    Rs = [0.01, 0.03, 0.10]
    NT = 24
    res = {}
    for R in Rs:
        for K in Ks:
            rng = np.random.default_rng(500)
            full, tri, ana = [], [], []
            for t in range(NT):
                f, tr, an = one_layout(K, R, rng)
                full.append(f); tri.append(tr); ana.append(an)
            res[(R, K)] = (np.concatenate(full), np.concatenate(tri), np.concatenate(ana))
            f, tr, an = res[(R, K)]
            print('R=%4.0f%% K=%2d  max|e_full| %6.2f  max T %6.2f  max B %8.1f  B/e %5.2f  violations T %5.1f%%  B %5.1f%%'
                  % (100 * R, K, f.max(), tr.max(), an.max(), an.max() / f.max(), 100 * np.mean(f > tr * 1.0001 + 0.05), 100 * np.mean(f > an + 0.05)))

    np.savez('out/s55_results.npz',
             err=np.concatenate([res[(R, K)][0] for R in Rs for K in Ks]),
             bound=np.concatenate([res[(R, K)][2] for R in Rs for K in Ks]),
             R=np.concatenate([np.full(res[(R, K)][0].size, R) for R in Rs for K in Ks]),
             K=np.concatenate([np.full(res[(R, K)][0].size, K) for R in Rs for K in Ks]))
    fig, ax = plt.subplots(1, 2, figsize=(3.5, 1.9), layout='constrained')
    cols = {0.01: FS.C_GOOD, 0.03: FS.GREEN, 0.10: FS.C_MEAS}
    mk = {0.01: 's', 0.03: '^', 0.10: 'o'}
    for R in Rs:
        f = np.concatenate([res[(R, K)][0] for K in Ks])
        tr = np.concatenate([res[(R, K)][2] for K in Ks])
        ax[0].loglog(tr, f, mk[R], ms=2.2, mfc='none', mec=cols[R], mew=0.6, alpha=0.7, label='$R=%g\\%%$' % (100 * R))
    lim = (0.05, 3000)
    ax[0].plot(lim, lim, '-', color=FS.C_THEORY, lw=0.9)
    ax[0].set_xlim(lim); ax[0].set_ylim(0.05, 300)
    ax[0].set_xlabel('Bound $\delta\lambda_k^{\max}$ [pm]')
    ax[0].set_ylabel('Full-model error [pm]')
    ax[0].legend(fontsize=6, loc='upper left', handlelength=1.2, borderaxespad=0.3)
    FS.letter(ax[0], 'a')
    for R in Rs:
        fmax = [res[(R, K)][0].max() for K in Ks]
        tmax = [res[(R, K)][1].max() for K in Ks]
        amax = [res[(R, K)][2].max() for K in Ks]
        ax[1].semilogy(Ks, fmax, mk[R] + '-', color=cols[R], ms=3.2, lw=1.2, label='full, $R=%g\\%%$' % (100 * R))
        ax[1].semilogy(Ks, amax, ls=(0, (5, 2)), color=cols[R], lw=1.0)
    ax[1].plot([], [], ls=(0, (5, 2)), color=FS.C_THEORY, lw=1.0, label='analytic bound')
    ax[1].set_xlabel('Gratings on the fiber, $K$')
    ax[1].set_ylabel('Largest error [pm]')
    ax[1].set_ylim(1, 2000)
    ax[1].legend(fontsize=5.5, loc='upper left', handlelength=1.8, borderaxespad=0.3, labelspacing=0.2)
    FS.letter(ax[1], 'b')
    fig.savefig('figs/fig_s55_error_bound.pdf'); fig.savefig('figs/fig_s55_error_bound.png', dpi=300)
    print('saved figs/fig_s55_error_bound.pdf and .png')
