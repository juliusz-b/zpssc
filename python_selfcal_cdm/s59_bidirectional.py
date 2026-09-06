"""s59_bidirectional.py - deshadowing from both fiber ends, in closed form.

Read from the source end, grating k shows S_k^f = R_k prod_{j<k} T_j^2.
Read from the far end (a loop, a switch, or a second circulator), it shows
S_k^b = R_k prod_{j>k} T_j^2, with T_j = 1 - R_j. Their product is
R_k^2 (T_tot / T_k)^2, where T_tot = prod_j T_j is the transmission of the
whole array, measured once. Hence

    R_k / (1 - R_k) = sqrt(S_k^f S_k^b) / T_tot,   pointwise in lambda,

which recovers every R_k(lambda) without the sequential division of (14) and
without its error accumulation. Only three measured quantities enter, and
none of them is a fitted line.

Test: K gratings sharing one band, random detunings +-200 pm, Gaussian lines
of 250 pm, relative noise on every measured spectrum. Compare the centre
error of (i) the raw forward reading, (ii) sequential deshadowing with fitted
transmissions and the 0.05 floor (s12/s15), (iii) the bidirectional closed
form. Spectral model only (no code, no chirp), so the shadowing mechanism is
isolated.

Output: out/s59_bidirectional.txt
"""
import os
import numpy as np
import common as C

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
W = 200.0 / PM
LAM = np.linspace(-650.0, 650.0, 64) / PM


def peel(S, floor=0.05):
    K = S.shape[0]; T = np.ones(LAM.size); cent = np.empty(K)
    for k in range(K):
        y = S[k] / np.maximum(T, floor)
        a, mu, sg, _ = C.gauss_fit_full(LAM, y)
        cent[k] = mu
        T = T * (1.0 - np.clip(a, 0, 0.99) * np.exp(-0.5 * ((LAM - mu) / max(sg, 1e-3)) ** 2)) ** 2
    return cent


def trial(K, R, noise_rel, seed):
    rng = np.random.default_rng(seed)
    nub = rng.uniform(-W, W, K)
    Rk = R * np.exp(-0.5 * ((LAM[None, :] - nub[:, None]) / SIG) ** 2)
    T = 1.0 - Rk
    cumf = np.cumprod(np.concatenate([np.ones((1, LAM.size)), T[:-1] ** 2], axis=0), axis=0)   # prod_{j<k} T_j^2
    cumb = np.cumprod(np.concatenate([np.ones((1, LAM.size)), T[:0:-1] ** 2], axis=0), axis=0)[::-1]  # prod_{j>k} T_j^2
    Sf = Rk * cumf; Sb = Rk * cumb; Ttot = np.prod(T, axis=0)
    nz = lambda X: X * (1.0 + rng.normal(0.0, noise_rel, X.shape)) + rng.normal(0.0, noise_rel * R, X.shape)
    Sf_m, Sb_m, Ttot_m = nz(Sf), nz(Sb), Ttot * (1.0 + rng.normal(0.0, noise_rel, Ttot.shape))
    q = np.sqrt(np.clip(Sf_m * Sb_m, 0, None)) / np.clip(Ttot_m, 1e-6, None)     # R/(1-R)
    R_bi = q / (1.0 + q)
    fit = lambda X: np.array([C.gauss_fit_peak(LAM, X[k]) for k in range(K)])
    e_raw = (fit(Sf_m) - nub) * PM
    e_seq = (peel(Sf_m) - nub) * PM
    e_bi = (fit(R_bi) - nub) * PM
    # amplitude recovery, relative RMS error of the peak reflectivity
    amp_err = np.sqrt(np.mean(((R_bi.max(axis=1) - Rk.max(axis=1)) / R) ** 2))
    return e_raw, e_seq, e_bi, amp_err


if __name__ == '__main__':
    lines = []
    say = lambda s='': (print(s), lines.append(s))
    rms = lambda e: float(np.sqrt(np.mean(e ** 2)))
    for R in (0.10, 0.03):
        for noise in (0.0, 1e-3, 1e-2):
            say('=== R = %.0f percent, relative noise %.0e, 12 trials' % (100 * R, noise))
            say('%6s %10s %12s %14s %14s' % ('K', 'raw', 'sequential', 'bidirectional', 'R error (bi)'))
            for K in (4, 8, 16, 32, 48, 64):
                out = [trial(K, R, noise, 100 + t) for t in range(12)]
                er = np.concatenate([o[0] for o in out]); es = np.concatenate([o[1] for o in out]); eb = np.concatenate([o[2] for o in out])
                ae = np.mean([o[3] for o in out])
                say('%6d %10.2f %12.2f %14.2f %13.1f%%' % (K, rms(er), rms(es), rms(eb), 100 * ae))
            say()
    os.makedirs('out', exist_ok=True)
    open('out/s59_bidirectional.txt', 'w').write('\n'.join(lines) + '\n')
    print('saved out/s59_bidirectional.txt')
