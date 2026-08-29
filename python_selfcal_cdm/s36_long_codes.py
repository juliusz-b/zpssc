"""s36_long_codes.py - leakage ceiling for the long codes, N = 1023 and 2047.

The sweep in s22 stops at K = 192, so N = 1023 printed as '>= 192' and the
invariant panel of Fig. 9 had nothing to show past N = 511. The ceiling is
supposed to grow like 0.37N, which puts N = 1023 near 380 and N = 2047 near
760, both well outside that ladder. This script runs the same leakage model
on a ladder that reaches them, so the two extra points on the invariant are
simulated rather than extrapolated from the very law the panel is testing.

Same model as s22: gratings on random delay bins, m-sequence side lobes,
detector noise, Gaussian fit, RMS over the array.
"""
import numpy as np
import common as C

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
DETUNE = 25.0
R = 0.10
M_GRID = 256
nu = np.linspace(-2.6 * F, 2.6 * F, M_GRID)
SIGMA_N = 1e-4
TARGET = 10.0

B = 25e6
FS_ADC = 20e6
OSR = 4
M_SLOTS = 64
C_LIGHT = 2.998e8
N_GROUP = 1.468


def acorr(nbits):
    m = 1.0 - 2.0 * C._mls01(nbits)
    return C.periodic_xcorr(m, m)


def leak_error(K, nchips, ac, seed=0, ntrials=3):
    out = []
    for t in range(ntrials):
        r = np.random.default_rng(seed + t)
        b = np.sort(r.choice(np.arange(1, nchips), size=K, replace=False))
        nub = r.uniform(-DETUNE, DETUNE, size=K)
        A = R * np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
        W = ac[(b[:, None] - b[None, :]) % nchips]
        np.fill_diagonal(W, 0.0)
        S = A + W @ A + r.normal(0, SIGMA_N, A.shape)
        e = [(C.gauss_fit_peak(nu, S[k]) - nub[k]) * PM for k in range(K)]
        out.append(np.sqrt(np.mean(np.array(e) ** 2)))
    return float(np.mean(out))


def k_at_target(nchips, ac, ladder):
    ks = np.array([k for k in ladder if k < nchips])
    e = np.array([leak_error(int(k), nchips, ac, seed=1000 + int(k))
                  for k in ks])
    for k, err in zip(ks, e):
        print('     K = %4d -> %7.2f pm' % (k, err))
    if e[-1] <= TARGET:
        return float(ks[-1]), True
    i = int(np.argmax(e > TARGET))
    if i == 0:
        return 0.0, False
    x0, x1, y0, y1 = ks[i - 1], ks[i], e[i - 1], e[i]
    return float(x0 + (TARGET - y0) * (x1 - x0) / (y1 - y0)), False


def ets_step(nchips):
    t_seq = nchips / B
    per = max(1, int(np.floor(FS_ADC * t_seq)))
    return int(np.ceil(nchips * OSR / per)) * t_seq


LADDER = [64, 128, 192, 256, 320, 384, 448, 512, 640, 768, 896]

print('leakage ceiling and refresh rate for the long codes')
print('(same model as s22, ladder extended to K = 896)')
for nbits in (10, 11):
    N = 2 ** nbits - 1
    print('  N = %d:' % N)
    k, sat = k_at_target(N, acorr(nbits), LADDER)
    frame = M_SLOTS * ets_step(N)
    fr = 1.0 / frame
    print('   -> K at %.0f pm: %.1f%s' % (TARGET, k, ' (ladder saturated)' if sat else ''))
    print('      frame %.2f ms, refresh %.1f Hz, K*f_r = %.0f sensor*Hz'
          % (frame * 1e3, fr, k * fr))
    print('      for reference, 0.37N = %.0f' % (0.37 * N))
