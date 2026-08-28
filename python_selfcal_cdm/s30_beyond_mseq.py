"""s30_beyond_mseq.py - is there a code between the m-sequence and the Golay pair?

Question from the design session of 2026-08-28: can an optimization produce a
sequence (or a small set, e.g. a triple) that beats the m-sequence's -1/N
periodic side lobe without paying the Golay pair's second transmission?

The stepped receiver integrates whole sequence periods, so the figure of merit
is the PERIODIC autocorrelation side lobe (Section II-B of the paper). Four
routes are checked.

1. PARITY BOUNDS (exact, no simulation needed). For a binary +-1 sequence of
   odd length N every periodic autocorrelation value R(tau) = N - 2d(tau) has
   the parity of N, so |R| >= 1: the m-sequence, which achieves exactly -1 at
   every lag, is optimal among single binary sequences. Modulo 4 the values
   satisfy R == N (mod 4). For N == 3 (mod 4), e.g. N = 127:
     - a PAIR sums to R1+R2 == 6 == 2 (mod 4), so |sum| >= 2 and the
       normalized floor 2/(2N) = 1/N: two shots buy nothing.
     - a TRIPLE sums to == 9 == 1 (mod 4), so |sum| >= 1 and the floor is
       1/(3N): three times lower leakage for three times the acquisition.
   The bit-flip search below shows how close a real triple gets.
   For comparison, the Golay pair (even N = 128) reaches exactly zero in two
   shots, so any odd-N triple is dominated: more time, more leakage.

2. TRIPLE OPTIMIZATION. Greedy bit-flip descent on max |sum of the three
   periodic autocorrelations|, ISL as tie-break, random restarts plus an
   m-sequence-seeded start.

3. RECEIVER-SIDE ALTERNATIVE: the cyclic inverse (mismatched) filter. The
   unipolar m-sequence record has a flat power spectrum with no zero bins
   (DC = (N+1)/2, elsewhere |C|^2 = (N+1)/4), so the cyclic filter
   r = IDFT(1/C*) exists and gives EXACTLY zero periodic side lobes with one
   transmission. Cost 1: a small noise penalty, measured below in dB against
   the paper's mean-removed matched receiver. Cost 2: the filter uses the DC
   bin, so a constant record offset is no longer rejected. Measured below: the
   offset leaks as a delay-independent pedestal, which the baseline term of
   the spectral fit removes anyway. POSTSCRIPT, same session: this filter is
   not new and has a name. For an m-sequence the cyclic inverse coincides,
   up to the scale (N+1)/2, with the bipolar replica in the bit-1 -> +1
   convention (replica sum = +1) used WITHOUT mean removal. That pairing is
   periodic non-coherent pulse compression: Ipatov 1992, Levanon et al.,
   IET RSN 10(1):216, 2016, Jahangir & Ali, Aerosp. Sci. Technol. 53:188,
   2016 (pointer: Bhatt, Trait. Signal 37(3):477, 2020, Sec. 3). The -1/N
   floor of the paper's receiver is created by the mean removal itself,
   given the correctly signed replica. Verified in section 3b below.

4. ANALOG DRIVE. Intensity modulation is not restricted to binary levels. An
   alternating-projection search (flat spectrum <-> box constraint 0..1) looks
   for a nonnegative waveform with perfect periodic autocorrelation. Reported
   with its modulation depth. Caveat for the paper: a multilevel drive feeds
   the adiabatic chirp term at every level, so this route trades a code error
   for a source error and is NOT recommended without a chirp budget.

Common yardstick: the s23-style leakage-only Bragg error at K = 32, N = 127,
R = 10%, 20 random layouts, same estimator as Table I of the paper.
"""
import numpy as np
from scipy.optimize import curve_fit
import warnings; warnings.filterwarnings('ignore')
import common as C

rng_global = np.random.default_rng(30)

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
W = 25.0
N = 127

# ---------------------------------------------------------------------------
# codes
# ---------------------------------------------------------------------------
ms01 = C._mls01(7)                       # unipolar {0,1}
ms = 1.0 - 2.0 * ms01                    # bipolar +-1 (0 -> +1, 1 -> -1)
ACORR_MS = C.periodic_xcorr(ms, ms)      # normalized, peak 1, side -1/N

gold = C.gold_set(7)
g3 = 1.0 - 2.0 * np.asarray(gold[3], dtype=float)
ACORR_GOLD = C.periodic_xcorr(g3, g3)

# Golay pair, length 128, standard doubling construction
a, b = np.array([1.0]), np.array([1.0])
for _ in range(7):
    a, b = np.concatenate([a, b]), np.concatenate([a, -b])
NG = len(a)
pr = lambda s: np.real(np.fft.ifft(np.abs(np.fft.fft(s)) ** 2))
GOLAY_SUM = (pr(a) + pr(b)) / (2.0 * NG)   # exactly delta


def peak_side(rho):
    return float(np.max(np.abs(rho[1:])))


# ---------------------------------------------------------------------------
# 1+2. binary search: single, and the triple
# ---------------------------------------------------------------------------
def descend(seqs, iters=200):
    """Greedy bit-flip descent on (max side lobe, ISL) of the summed
    periodic autocorrelation of the given bipolar sequences."""
    seqs = [s.copy() for s in seqs]
    S = [np.fft.fft(s) for s in seqs]
    om = [np.exp(-2j * np.pi * np.arange(N) * i / N) for i in range(N)]

    def metrics():
        P = np.sum([np.abs(Si) ** 2 for Si in S], axis=0)
        r = np.real(np.fft.ifft(P))
        side = np.abs(r[1:])
        return float(side.max()), float(np.sum(side ** 2))

    best = metrics()
    for _ in range(iters):
        improved = False
        for j in range(len(seqs)):
            for i in range(N):
                old = S[j].copy()
                S[j] = S[j] - 2.0 * seqs[j][i] * om[i]
                m = metrics()
                if m < best:
                    seqs[j][i] = -seqs[j][i]
                    best = m
                    improved = True
                else:
                    S[j] = old
        if not improved:
            break
    return best, seqs


print('=== 1. parity bounds at N = %d (values mod 4: R == 3) ===' % N)
print('  single: |R| >= 1  -> best possible 1/N = %.4f (m-sequence: %.4f)'
      % (1.0 / N, peak_side(np.abs(ACORR_MS))))
print('  pair  : |R1+R2| >= 2 -> floor 2/2N = 1/N: no gain from two shots')
print('  triple: |sum| >= 1 -> floor 1/3N = %.4f' % (1.0 / (3 * N)))
print()

best_single = min(descend([1.0 - 2.0 * (rng_global.random(N) < 0.5)])[0][0]
                  for _ in range(8)) / N
print('=== 2. bit-flip search ===')
print('  single, 8 restarts: best max side lobe %.4f (m-seq: %.4f) '
      '-> m-sequence not beaten, as parity demands' % (best_single, 1.0 / N))

cands = []
for t in range(12):
    seed = [1.0 - 2.0 * (rng_global.random(N) < 0.5) for _ in range(3)]
    cands.append(descend(seed))
cands.append(descend([ms.copy(), np.roll(ms, 43).copy(),
                      1.0 - 2.0 * (rng_global.random(N) < 0.5)]))
(best_m, best_isl), best_seqs = min(cands, key=lambda x: x[0])
TRIPLE_SUM = np.real(np.fft.ifft(np.sum(
    [np.abs(np.fft.fft(s)) ** 2 for s in best_seqs], axis=0))) / (3.0 * N)
print('  triple, 13 restarts: best max |sum|/3N = %.4f '
      '(parity floor %.4f, m-seq single %.4f)'
      % (best_m / (3 * N), 1.0 / (3 * N), 1.0 / N))
print('  trivial triple (m-seq three times): sum = -3 flat -> %.4f, '
      'no better than one shot' % (3.0 / (3 * N)))
print('  golay pair, N = 128, 2 shots: max side lobe %.1e (exact zero)'
      % peak_side(GOLAY_SUM))
print()

# ---------------------------------------------------------------------------
# 3. cyclic inverse filter on the unipolar m-sequence record
# ---------------------------------------------------------------------------
c01 = ms01.astype(float)                   # what the laser actually sends
Cf = np.fft.fft(c01)
assert np.min(np.abs(Cf)) > 1e-9, 'spectrum has a zero bin, no inverse'
r_inv = np.real(np.fft.ifft(1.0 / np.conj(Cf)))
RHO_INV = np.real(np.fft.ifft(np.fft.fft(c01) * np.fft.fft(r_inv).conj()))
# matched, mean-removed receiver of the paper (Section II-A)
c0 = c01 - c01.mean()
msb = 2.0 * ms01 - 1.0                     # bipolar replica, chip 1 -> +1
RHO_MATCH = np.real(np.fft.ifft(np.fft.fft(c0) * np.fft.fft(msb).conj()))
RHO_MATCH /= RHO_MATCH[0]

print('=== 3. cyclic inverse filter (one transmission) ===')
print('  side lobes: matched mean-removed %.4f, inverse %.1e (exact zero)'
      % (peak_side(RHO_MATCH), peak_side(np.abs(RHO_INV) / RHO_INV[0])))

# noise penalty, measured: unit echo at delay 40 + white noise
sig_n, trials = 0.05, 400
rng = np.random.default_rng(7)
y0 = np.roll(c01, 40)
corr_m = lambda y: np.real(np.fft.ifft(
    np.fft.fft(y - y.mean()) * np.fft.fft(msb).conj()))
corr_i = lambda y: np.real(np.fft.ifft(np.fft.fft(y) * np.fft.fft(r_inv).conj()))
zm0, zi0 = corr_m(y0), corr_i(y0)
nm, ni = [], []
for _ in range(trials):
    y = y0 + rng.normal(0.0, sig_n, N)
    nm.append(corr_m(y) - zm0); ni.append(corr_i(y) - zi0)
snr_m = zm0[40] / np.std(nm)
snr_i = zi0[40] / np.std(ni)
print('  noise penalty of the inverse: %.2f dB (20log10 convention)'
      % (20 * np.log10(snr_m / snr_i)))

# offset robustness: constant offset -> flat pedestal only
z_off = corr_i(y0 + 0.5)
ped = z_off - zi0
print('  0.5 offset -> pedestal mean %.4f, flatness (std) %.1e '
      '-> removed by the baseline term of the fit' % (ped.mean(), ped.std()))
print()

# ---------------------------------------------------------------------------
# 3b. the inverse filter has a name: the NCPC replica (Ipatov/Levanon)
# ---------------------------------------------------------------------------
xc_ncpc = np.real(np.fft.ifft(np.fft.fft(c01) * np.fft.fft(msb).conj()))
xm_mr = np.real(np.fft.ifft(np.fft.fft(c01 - c01.mean()) * np.fft.fft(msb).conj()))
print("=== 3b. equivalence: inverse filter == NCPC bipolar replica ===")
print('  replica sum = %+d (bit 1 -> +1), no mean removal:' % int(msb.sum()))
print('    peak %.0f = (N+1)/2, max side lobe %.1e (exact zero)'
      % (xc_ncpc[0], np.abs(xc_ncpc[1:]).max()))
print('  b vs (N+1)/2 * r_inv: max deviation %.1e (identical up to scale)'
      % np.abs(msb - (N + 1) / 2.0 * r_inv).max())
print('  same replica AFTER mean removal: side lobe %.6f, i.e. the -1/N'
      % (xm_mr[1] / xm_mr[0]))
print('  floor of the paper is created by the mean removal itself')
print()

# ---------------------------------------------------------------------------
# 4. analog nonnegative drive with perfect periodic autocorrelation
# ---------------------------------------------------------------------------
x = rng_global.random(N)
for _ in range(4000):
    X = np.fft.fft(x)
    mag = np.abs(X[1:]).mean()
    X[1:] = mag * X[1:] / np.maximum(np.abs(X[1:]), 1e-12)
    x = np.clip(np.real(np.fft.ifft(X)), 0.0, 1.0)
Xf = np.fft.fft(x)
rho_a = np.real(np.fft.ifft(np.abs(Xf) ** 2))
rho_a_side = peak_side(rho_a - np.mean(x) ** 2 * N) / (rho_a[0] - np.mean(x) ** 2 * N)
print('=== 4. analog drive (alternating projections, 0..1 box) ===')
print('  best AC-coupled side lobe %.1e, modulation depth min/max = '
      '%.2f/%.2f, mean level %.2f' % (rho_a_side, x.min(), x.max(), x.mean()))
print('  caveat: multilevel drive feeds adiabatic chirp at every level')
print()

# ---------------------------------------------------------------------------
# common yardstick: leakage-only Bragg error at K = 32 (s23-style trial)
# ---------------------------------------------------------------------------
K, R_B, M = 32, 0.10, 256
nu = np.linspace(-2.6 * F, 2.6 * F, M)


def lsq_full(xg, y, mu0):
    f = lambda t, aa, m, s, cc: aa * np.exp(-0.5 * ((t - m) / s) ** 2) + cc
    try:
        p, _ = curve_fit(f, xg, y, p0=[max(y.max(), 1e-9), mu0, SIG, 0.0],
                         maxfev=20000)
        return float(p[1])
    except Exception:
        return mu0


def leak_rms(rho, nch, trials=20):
    errs = []
    for s_ in range(trials):
        rng = np.random.default_rng(100 + s_)
        bins = np.sort(rng.choice(np.arange(1, nch), size=K, replace=False))
        nub = rng.uniform(-W, W, size=K)
        A = R_B * np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
        Wl = rho[(bins[:, None] - bins[None, :]) % nch]
        np.fill_diagonal(Wl, 0.0)
        S = A + Wl @ A
        mus = np.array([lsq_full(nu, S[k], nub[k]) for k in range(K)])
        errs.append(np.sqrt(np.mean(((mus - nub) * PM) ** 2)))
    return float(np.mean(errs))


print('=== leakage-only Bragg error, K = 32, R = 10%%, 20 layouts ===')
rows = [
    ('m-sequence, matched (paper)', RHO_MATCH, N, 1),
    ('Gold, matched', ACORR_GOLD, N, 1),
    ('optimized triple', TRIPLE_SUM, N, 3),
    ('Golay pair, N=128', GOLAY_SUM, NG, 2),
    ('m-sequence + inverse filter', np.where(np.arange(N) == 0, 1.0, 0.0), N, 1),
]
print('   %-30s %10s %6s' % ('scheme', 'RMS [pm]', 'shots'))
for name, rho, nch, shots in rows:
    print('   %-30s %10.2f %6d' % (name, leak_rms(rho, nch), shots))
print()
print('VERDICT: single binary cannot beat the m-sequence (parity). Pairs at')
print('odd N buy nothing, triples reach at best 1/3N and are dominated by the')
print('Golay pair. The real free lunch is receiver-side and has a name:')
print('periodic NCPC (Ipatov/Levanon), i.e. the correctly signed bipolar')
print('replica without mean removal. Zero periodic leakage, one transmission,')
print('no SNR cost, and no extra computation over the matched correlator.')
