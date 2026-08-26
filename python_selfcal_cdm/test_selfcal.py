"""test_selfcal.py - checks for the physics and the numerics this study rests on.

Run with:  python test_selfcal.py
No test framework is required; the file prints one line per check and exits with
a non-zero status if anything fails.

The checks are grouped by what they protect:

  CODES        the correlation bounds that set the despreading floor. If these
               break, every leakage number in s12 and s14 is wrong.
  ESTIMATORS   the peak fitters, including the full-parameter fit that the
               deshadowing recursion depends on.
  PHYSICS      the limiting behaviours the method claims: a symmetric chirp on a
               symmetric grating shifts nothing, an asymmetric one does, and a
               ghost built from identical gratings does not move the peak.
  ARRAY MODEL  the delay algebra behind the uniform-versus-randomised spacing
               result, and the monotonic growth of the error with array size.
  ACQUISITION  the spatial-resolution formula and the equivalent-time budget.
"""
import sys
import numpy as np
import warnings; warnings.filterwarnings('ignore')
import common as C

FAILED = []


def check(name, cond, detail=''):
    ok = bool(cond)
    print('  %-58s %s%s' % (name, 'PASS' if ok else 'FAIL',
                            ('   ' + detail) if detail else ''))
    if not ok:
        FAILED.append(name)


PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ

# ---------------------------------------------------------------------------
print('CODES')
# ---------------------------------------------------------------------------
g7 = C.gold_set(7)
auto7, cross7 = C.code_corr_stats(g7)
check('Gold n=7 cross-correlation meets the 17/127 bound',
      abs(cross7 - 17 / 127) < 1e-9, 'got %.6f' % cross7)
check('Gold n=7 autocorrelation side lobe equals the same bound',
      abs(auto7 - 17 / 127) < 1e-9, 'got %.6f' % auto7)

k6 = C.kasami_small(6)
a6, x6 = C.code_corr_stats(k6)
check('Kasami small n=6 cross-correlation meets the 9/63 bound',
      abs(x6 - 9 / 63) < 1e-9, 'got %.6f' % x6)

m7 = 1.0 - 2.0 * C._mls01(7)
ac = C.periodic_xcorr(m7, m7)
check('m-sequence autocorrelation side lobe is exactly -1/N',
      np.allclose(ac[1:], -1.0 / 127, atol=1e-9),
      'max deviation %.2e' % np.max(np.abs(ac[1:] + 1 / 127)))
check('m-sequence beats Gold on the side lobe by 17x',
      abs((17 / 127) / np.max(np.abs(ac[1:])) - 17.0) < 0.1)
check('m-sequence is balanced (DC term is one chip)',
      abs(m7.sum()) == 1.0, 'sum = %.1f' % m7.sum())


def golay_pair(k):
    a, b = np.array([1.0]), np.array([1.0])
    for _ in range(k):
        a, b = np.concatenate([a, b]), np.concatenate([a, -b])
    return a, b


ga, gb = golay_pair(7)
gsum = C.periodic_xcorr(ga, ga) + C.periodic_xcorr(gb, gb)
check('Golay pair: summed autocorrelation has no side lobes',
      np.allclose(gsum[1:], 0.0, atol=1e-10),
      'max side lobe %.2e' % np.max(np.abs(gsum[1:])))

# ---------------------------------------------------------------------------
print('ESTIMATORS')
# ---------------------------------------------------------------------------
x = np.linspace(-200.0, 200.0, 1200)
truth = 17.3
y = 0.42 * np.exp(-0.5 * ((x - truth) / 12.0) ** 2) + 0.05
check('gauss_fit_peak recovers a clean Gaussian centre to 0.01 GHz',
      abs(C.gauss_fit_peak(x, y) - truth) < 0.01,
      'got %.4f' % C.gauss_fit_peak(x, y))
a, mu, sg, c0 = C.gauss_fit_full(x, y)
check('gauss_fit_full recovers amplitude, centre, width and baseline',
      abs(a - 0.42) < 0.01 and abs(mu - truth) < 0.01 and abs(sg - 12.0) < 0.1
      and abs(c0 - 0.05) < 0.01,
      'a=%.3f mu=%.2f sig=%.2f c=%.3f' % (a, mu, sg, c0))
check('centroid and Gaussian fit agree on a symmetric line',
      abs(C.centroid(x, y, thr=0.3) - C.gauss_fit_peak(x, y)) < 0.5)

# ---------------------------------------------------------------------------
print('PHYSICS')
# ---------------------------------------------------------------------------
g = np.linspace(-6 * F, 6 * F, 1400)
sym = C.fmam_readout(g, 0.0, F, 0.4 * F, mean_off=0.0, skew=0.0, shape='tanh', asym=0.0)
check('symmetric chirp on a symmetric grating shifts nothing',
      abs(C.gauss_fit_peak(g, sym) * PM) < 0.2,
      '%.3f pm' % (C.gauss_fit_peak(g, sym) * PM))

skewed = C.fmam_readout(g, 0.0, F, 0.4 * F, mean_off=0.0, skew=1.2, shape='tanh',
                        asym=0.30)
base = C.fbg_tanh(g, 0.0, F, n_side=0.30)
shift = (C.gauss_fit_peak(g, skewed) - C.gauss_fit_peak(g, base)) * PM
check('skewed chirp on an asymmetric grating does shift the peak',
      abs(shift) > 1.0, '%.2f pm' % shift)

shifts = []
for ratio in (0.1, 0.2, 0.4):
    s_ = C.fmam_readout(g, 0.0, F, ratio * F, mean_off=0.0, skew=1.2, shape='tanh',
                        asym=0.30)
    shifts.append(abs((C.gauss_fit_peak(g, s_) - C.gauss_fit_peak(g, base)) * PM))
check('chirp-induced shift grows with the chirp excursion',
      shifts[0] < shifts[1] < shifts[2],
      '%.2f < %.2f < %.2f pm' % tuple(shifts))

pure_off = C.fmam_readout(g, 0.0, F, 1e-6, mean_off=2.0, skew=0.0, shape='gauss')
check('a pure frequency offset moves the peak by exactly minus that offset',
      abs(C.gauss_fit_peak(g, pure_off) + 2.0) < 0.02,
      '%.4f GHz' % C.gauss_fit_peak(g, pure_off))

prim = C.fbg_tanh(g, 0.0, F, 0.0)
for alpha in (0.05, 0.2):
    comb = prim + alpha * prim
    check('ghost from identical gratings does not move the peak (alpha=%.2f)' % alpha,
          abs(C.gauss_fit_peak(g, comb) * PM) < 0.05,
          '%.4f pm' % (C.gauss_fit_peak(g, comb) * PM))

# ---------------------------------------------------------------------------
print('ARRAY MODEL')
# ---------------------------------------------------------------------------
K, D, tau0 = 12, 7, 3
bins_u = tau0 + D * np.arange(K)


def ghost_stats(bins):
    """Third-order ghosts that fall inside the span occupied by the array, and
    how many of them coincide with a real delay bin."""
    inside = hit = 0
    lo, hi = bins.min(), bins.max()
    for ai in range(len(bins)):
        for bi in range(len(bins)):
            for ci in range(len(bins)):
                if bi < ai and bi < ci:
                    gb = bins[ai] - bins[bi] + bins[ci]
                    if lo <= gb <= hi:
                        inside += 1
                        hit += int(gb in bins)
    return inside, hit


ins_u, hit_u = ghost_stats(bins_u)
check('uniform spacing: every ghost inside the array span hits a real bin',
      hit_u == ins_u and ins_u > 0, '%d of %d' % (hit_u, ins_u))

fracs = []
for seed in range(8):
    rng = np.random.default_rng(seed)
    br = np.sort(rng.choice(np.arange(1, 127), size=K, replace=False))
    ins_r, hit_r = ghost_stats(br)
    fracs.append(hit_r / max(ins_r, 1))
frac_r = float(np.mean(fracs))
check('randomised spacing: only a small fraction of them do',
      frac_r < 0.25 and frac_r < 0.25 * (hit_u / ins_u),
      'uniform 100%%, randomised %.0f%% (mean of 8 layouts, worst %.0f%%)'
      % (100 * frac_r, 100 * max(fracs)))

SIG = F / 2.35482
nu = np.linspace(-2.6 * F, 2.6 * F, 64)


def leak_error(Kg, nchips, acorr, seed=0):
    r = np.random.default_rng(seed)
    b = np.sort(r.choice(np.arange(1, nchips), size=Kg, replace=False))
    nub = r.uniform(-25, 25, size=Kg)
    A = 0.1 * np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
    W = acorr[(b[:, None] - b[None, :]) % nchips]; np.fill_diagonal(W, 0.0)
    S = A + W @ A
    e = np.array([(C.gauss_fit_peak(nu, S[k]) - nub[k]) * PM for k in range(Kg)])
    return float(np.sqrt(np.mean(e ** 2)))


e8 = np.mean([leak_error(8, 127, ac, s) for s in range(4)])
e32 = np.mean([leak_error(32, 127, ac, s) for s in range(4)])
check('code leakage error grows with array size',
      e32 > 2 * e8, 'K=8: %.2f pm, K=32: %.2f pm' % (e8, e32))

m9 = 1.0 - 2.0 * C._mls01(9); ac9 = C.periodic_xcorr(m9, m9)
e32_long = np.mean([leak_error(32, 511, ac9, s) for s in range(4)])
check('a four times longer code cuts the leakage error about fourfold',
      2.5 < e32 / e32_long < 6.0,
      'N=127: %.2f pm, N=511: %.2f pm, ratio %.2f' % (e32, e32_long, e32 / e32_long))

gold_ac = C.periodic_xcorr(g7[3], g7[3])
e32_gold = np.mean([leak_error(32, 127, gold_ac, s) for s in range(4)])
check('m-sequence beats Gold in the swept (single-code) architecture',
      e32_gold > 2 * e32, 'Gold %.2f pm vs m-seq %.2f pm' % (e32_gold, e32))

shp = np.exp(-0.5 * ((nu[None, :] - np.array([[-15.0], [12.0]])) / SIG) ** 2)
t1 = (1.0 - 0.10 * shp[0]) ** 2
check('spectral shadowing biases a downstream grating',
      abs((C.gauss_fit_peak(nu, 0.1 * shp[1] * t1) - 12.0) * PM) > 1.0,
      '%.2f pm' % ((C.gauss_fit_peak(nu, 0.1 * shp[1] * t1) - 12.0) * PM))
check('shadowing bias vanishes when the upstream grating is far detuned',
      abs((C.gauss_fit_peak(nu, 0.1 * shp[1] *
                            (1 - 0.10 * np.exp(-0.5 * ((nu + 200) / SIG) ** 2)) ** 2)
           - 12.0) * PM) < 0.2)

# ---------------------------------------------------------------------------
print('BANDING AND COLLISIONS')
# ---------------------------------------------------------------------------
BAND_STEP = 8.0 * F
sig_ = F / 2.35482
overlap = np.exp(-0.5 * (BAND_STEP / sig_) ** 2)
check('gratings one WDM band apart have negligible spectral overlap',
      overlap < 1e-12, 'exp term %.1e at 2 nm pitch' % overlap)

rngc = np.random.default_rng(5)
Kc, Nc, TRIALS = 32, 511, 4000
hits = 0
for _ in range(TRIALS):
    b = rngc.integers(0, Nc, size=Kc)
    hits += int(len(np.unique(b)) < Kc)
p_mc = hits / TRIALS
i = np.arange(1, Kc)
p_th = 1.0 - np.prod(1.0 - i / Nc)
check('birthday formula matches Monte Carlo for co-binning probability',
      abs(p_mc - p_th) < 0.03, 'theory %.3f, Monte Carlo %.3f' % (p_th, p_mc))

pairs_mc = np.mean([np.sum(np.bincount(rngc.integers(0, Nc, size=Kc),
                                       minlength=Nc) >= 2) for _ in range(2000)])
pairs_th = Kc * (Kc - 1) / (2.0 * Nc)
check('expected co-binned pairs close to K(K-1)/2N',
      abs(pairs_mc - pairs_th) < 0.12,
      'theory %.2f, Monte Carlo %.2f' % (pairs_th, pairs_mc))

# ---------------------------------------------------------------------------
print('ACQUISITION')
# ---------------------------------------------------------------------------
c_light, n_g = 2.99792458e8, 1.468
dz = lambda B: c_light / (2 * n_g * B)
check('25 Mchip/s gives about 4.1 m of spatial resolution',
      abs(dz(25e6) - 4.086) < 0.02, '%.3f m' % dz(25e6))
check('doubling the chip rate halves the resolution cell',
      abs(dz(50e6) / dz(25e6) - 0.5) < 1e-9)


def ets(B, fs, nchips=127, osr=4):
    t_seq = nchips / B
    per = max(1, int(np.floor(fs * t_seq)))
    return int(np.ceil(nchips * osr / per)) * t_seq


check('equivalent-time acquisition at 25 Mchip/s with a 3.6 MS/s converter'
      ' stays under 1 ms', ets(25e6, 3.6e6) < 1e-3, '%.1f us' % (ets(25e6, 3.6e6) * 1e6))
check('a faster converter shortens the equivalent-time acquisition',
      ets(25e6, 20e6) < ets(25e6, 3.6e6),
      '%.1f us vs %.1f us' % (ets(25e6, 20e6) * 1e6, ets(25e6, 3.6e6) * 1e6))

# ---------------------------------------------------------------------------
print('UNITS')
# ---------------------------------------------------------------------------
check('1 GHz corresponds to about 8 pm at 1550 nm', abs(PM - 8.0) < 0.001,
      '%.4f pm/GHz' % PM)
check('the procured 250 pm linewidth is 31.25 GHz',
      abs(F - 31.25) < 1e-6, '%.4f GHz' % F)
check('the Peltier stability translates to about 1 pm',
      abs(C.TEMP_STAB_C * C.TEMP_COEF_PM_PER_C - 1.0) < 1e-9)

print()
if FAILED:
    print('%d CHECK(S) FAILED: %s' % (len(FAILED), '; '.join(FAILED)))
    sys.exit(1)
print('all checks passed')
