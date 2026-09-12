"""Deshadowing in the power model: error of the LAST grating after sequential correction.

(a) RMS error of the last grating (over random layouts) against the number of gratings K, for several
    peak reflectivities R, with and without deshadowing. Same model as s12_capacity (shadowing, ghosts,
    m-sequence leakage, detector noise, local Gaussian fit), randomized delay bins, N = 127.
(b) Robustness to a spread of reflectivities: every grating gets R_k = R (1 + u), u uniform in
    [-s, s], and the correction uses either the fitted line amplitude or the nominal R. RMS error of
    the last grating against the spread s, for K = 16.

Output: figs/fig_s56_peel_depth.{pdf,png}
"""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS
FS.apply()

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
NBITS, N_CHIPS, M_SLOTS = 7, 127, 64
BAND_GHZ, DETUNE_GHZ = 2.6 * F, 25.0
nu = np.linspace(-BAND_GHZ, BAND_GHZ, M_SLOTS)
_MSEQ = 1.0 - 2.0 * C._mls01(NBITS)
ACORR = C.periodic_xcorr(_MSEQ, _MSEQ)
P_SOURCE_W, LOSS_DB, NEP_W_RTHZ, CHIP_RATE_HZ, UNIPOLAR = 1.0e-3, 4.0, 0.5e-12, 25.0e6, 0.5
SIGMA_N = (NEP_W_RTHZ * np.sqrt(CHIP_RATE_HZ) / (P_SOURCE_W * 10 ** (-LOSS_DB / 10))) / (np.sqrt(N_CHIPS) * UNIPOLAR)


def layout(K, rng):
    return np.sort(rng.choice(np.arange(1, N_CHIPS), size=K, replace=False))


def spectra(K, Rk, rng):
    nub = rng.uniform(-DETUNE_GHZ, DETUNE_GHZ, size=K)
    shapes = np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
    tcum = np.ones((K, M_SLOTS))
    trans = (1.0 - Rk[:, None] * shapes) ** 2
    for k in range(1, K):
        tcum[k] = tcum[k - 1] * trans[k - 1]
    return nub, shapes, tcum


def ghosts(K, nub, Rk, bins, tcum):
    idx = np.arange(K)
    a, b, c = np.meshgrid(idx, idx, idx, indexing='ij')
    m = (b < a) & (b < c)
    a, b, c = a[m], b[m], c[m]
    gbin = bins[a] - bins[b] + bins[c]
    pos = np.clip(np.searchsorted(bins, gbin), 0, K - 1)
    hit = bins[pos] == gbin
    out = np.zeros((K, M_SLOTS))
    if not hit.any():
        return out
    a, b, c, target = a[hit], b[hit], c[hit], pos[hit]
    cbar = (nub[a] + nub[b] + nub[c]) / 3.0
    spread = (nub[a] - cbar) ** 2 + (nub[b] - cbar) ** 2 + (nub[c] - cbar) ** 2
    amp = Rk[a] * Rk[b] * Rk[c] * np.exp(-0.5 * spread / SIG ** 2)
    up = np.minimum(np.minimum(a, b), c)
    for i in range(len(amp)):
        out[target[i]] += amp[i] * tcum[up[i]] * np.exp(-0.5 * ((nu - cbar[i]) / (SIG / np.sqrt(3.0))) ** 2)
    return out


def measured(K, Rk, rng):
    bins = layout(K, rng)
    nub, shapes, tcum = spectra(K, Rk, rng)
    W = ACORR[(bins[:, None] - bins[None, :]) % N_CHIPS]; np.fill_diagonal(W, 0.0)
    direct = Rk[:, None] * shapes * tcum
    S = direct + W @ direct + ghosts(K, nub, Rk, bins, tcum)
    return nub, S + rng.normal(0, SIGMA_N, size=S.shape)


def peel(S, R0=None, floor=0.05):
    """Sequential deshadowing with the fitted amplitude (R0 None) or the nominal reflectivities."""
    K = S.shape[0]
    T = np.ones(M_SLOTS); cent = np.empty(K)
    for k in range(K):
        a, mu, sg, _ = C.gauss_fit_full(nu, S[k] / np.maximum(T, floor))
        cent[k] = mu
        amp = np.clip(a, 0.0, 0.99) if R0 is None else R0[k]
        T = T * (1.0 - amp * np.exp(-0.5 * ((nu - mu) / max(sg, 1e-3)) ** 2)) ** 2
    return cent


def last_error(K, R, spread, mode, ntrials, seed):
    """RMS error [pm] of the last grating over ntrials layouts. mode: 'raw', 'peel_fit', 'peel_R0'."""
    rng = np.random.default_rng(seed)
    errs = []
    for _ in range(ntrials):
        Rk = R * (1.0 + rng.uniform(-spread, spread, size=K))
        nub, S = measured(K, Rk, rng)
        if mode == 'raw':
            c = C.gauss_fit_peak(nu, S[-1])
        elif mode == 'peel_fit':
            c = peel(S)[-1]
        else:
            c = peel(S, R0=np.full(K, R))[-1]
        errs.append((c - nub[-1]) * PM)
    return float(np.sqrt(np.mean(np.square(errs))))


if __name__ == '__main__':
    Ks = [4, 8, 16, 24, 32, 48]
    Rs = [0.01, 0.03, 0.10, 0.20]
    NT = 24
    fig, ax = plt.subplots(1, 2, figsize=(3.5, 1.9), layout='constrained')
    cols = {0.01: FS.C_GOOD, 0.03: FS.GREEN, 0.10: FS.C_MEAS, 0.20: FS.C_GHOST}
    mk = {0.01: 's', 0.03: '^', 0.10: 'o', 0.20: 'D'}
    for R in Rs:
        raw = [last_error(K, R, 0.0, 'raw', NT, 100 + K) for K in Ks]
        pl = [last_error(K, R, 0.0, 'peel_R0', NT, 100 + K) for K in Ks]
        print('R=%3.0f%%  raw  %s' % (100 * R, np.round(raw, 1)))
        print('        peel %s' % np.round(pl, 1))
        ax[0].semilogy(Ks, raw, mk[R] + ':', color=cols[R], ms=3, lw=1.0, mfc='none')
        ax[0].semilogy(Ks, pl, mk[R] + '-', color=cols[R], ms=3, lw=1.3, label='$R=%g\\%%$' % (100 * R))
    ax[0].plot([], [], ':', color=FS.C_THEORY, lw=1.0, label='raw'); ax[0].plot([], [], '-', color=FS.C_THEORY, lw=1.3, label='deshadowed')
    ax[0].axhline(10, color='.5', ls='--', lw=.7)
    ax[0].set(xlabel='Gratings on the fiber, $K$', ylabel='RMS error, last grating [pm]', ylim=(0.5, 3000))
    ax[0].legend(fontsize=5, loc='upper left', handlelength=2.4, borderaxespad=0.3, labelspacing=0.15, ncol=2, columnspacing=0.8)
    FS.letter(ax[0], 'a')
    spreads = [0.0, 0.1, 0.2, 0.3, 0.5]
    K = 16
    for R in (0.03, 0.10):
        e_fit = [last_error(K, R, s, 'peel_fit', NT, 300) for s in spreads]
        e_R0 = [last_error(K, R, s, 'peel_R0', NT, 300) for s in spreads]
        e_raw = [last_error(K, R, s, 'raw', NT, 300) for s in spreads]
        print('R=%3.0f%% spread %s: raw %s | peel(R0) %s | peel(fit) %s' % (100 * R, spreads, np.round(e_raw, 1), np.round(e_R0, 1), np.round(e_fit, 1)))
        ax[1].plot(np.array(spreads) * 100, e_R0, mk[R] + '-', color=cols[R], ms=3, lw=1.3, label='$R=%g\\%%$, nominal $R_0$' % (100 * R))
        ax[1].plot(np.array(spreads) * 100, e_fit, mk[R] + '--', color=cols[R], ms=3, lw=1.0, mfc='none', label='$R=%g\\%%$, fitted $R_0$' % (100 * R))
    ax[1].set(xlabel='Reflectivity spread [%]', ylabel='RMS error, last of 16 [pm]', ylim=(4, 36))
    ax[1].legend(fontsize=5, loc='upper left', handlelength=2.4, borderaxespad=0.3, labelspacing=0.15)
    FS.letter(ax[1], 'b')
    fig.savefig('figs/fig_s56_peel_depth.pdf'); fig.savefig('figs/fig_s56_peel_depth.png', dpi=300)
    print('saved figs/fig_s56_peel_depth.pdf and .png')
