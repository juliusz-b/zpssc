"""s14_code_families.py - which spreading sequence should a swept CDM
interrogator use?

The reflex answer is Gold, because Gold sets are large and have a bounded
cross-correlation. That answer is right for a radio system in which many users
transmit at once. It is the wrong answer here, and the reason is the sweep.

With a swept source only one wavelength, hence one code, is on the fiber inside a
correlation window. The gratings are separated from each other by the code
AUTOCORRELATION side lobes, not by cross-correlation between codes. Maximal
length sequences have a periodic autocorrelation side lobe of exactly -1/N. Gold
sequences pay for their large family with three-valued autocorrelation reaching
(2^((n+1)/2)+1)/N, which is 17/127 for n = 7, i.e. 17 times worse. Golay
complementary pairs cancel their side lobes exactly, but only when the two
transmissions are added, so they cost twice the acquisition time.

What each property buys:
  * autocorrelation side lobe  -> grating-to-grating leakage inside one sweep
    slot. This is the term that matters for a swept source.
  * cross-correlation          -> only matters if several wavelengths are on the
    fiber simultaneously, e.g. a multi-emitter or comb source.
  * family size                -> number of wavelength slots addressable in
    parallel. For a swept source one code per slot in sequence is enough, so a
    family of 18 m-sequences is not a limitation.

Outputs
  (a) periodic autocorrelation of the candidates near zero lag.
  (b) RMS Bragg error vs array size for each family, shadowing and ghosts off so
      that only the code contributes.
"""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS
FS.apply()

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
M_SLOTS = 64
BAND_GHZ = 2.6 * F
DETUNE_GHZ = 25.0
R = 0.10
SIGMA_N = 1.1e-6            # same link-budget noise floor as s12

nu = np.linspace(-BAND_GHZ, BAND_GHZ, M_SLOTS)


def mseq_family(n):
    """All distinct m-sequences of degree n, built by decimation of one of them."""
    base = C._mls01(n)
    N = 2 ** n - 1
    seen, out = set(), []
    for q in range(1, N):
        if np.gcd(q, N) != 1:
            continue
        s = C._decimate01(base, q)
        key = min(tuple(np.roll(s, k)) for k in range(N))
        if key in seen:
            continue
        seen.add(key)
        out.append(1.0 - 2.0 * s)
    return np.array(out)


def golay_pair(k):
    """Complementary pair of length 2^k."""
    a, b = np.array([1.0]), np.array([1.0])
    for _ in range(k):
        a, b = np.concatenate([a, b]), np.concatenate([a, -b])
    return a, b


def acorr(c):
    return C.periodic_xcorr(c, c)


# --- candidate families -----------------------------------------------------
mseq7 = mseq_family(7)
gold7 = C.gold_set(7)
kas6 = C.kasami_small(6)
ga, gb = golay_pair(7)                       # length 128

FAMILIES = []
FAMILIES.append(('m-sequence, N=127', 127, len(mseq7), acorr(mseq7[0])))
FAMILIES.append(('Gold, N=127', 127, len(gold7), acorr(gold7[3])))
FAMILIES.append(('Kasami small, N=63', 63, len(kas6), acorr(kas6[1])))
FAMILIES.append(('Golay pair, N=128 (2 shots)', 128, 2, (acorr(ga) + acorr(gb)) / 2.0))


def cross_max(codes, limit=24):
    codes = codes[:limit]
    m = 0.0
    for i in range(len(codes)):
        for j in range(i + 1, len(codes)):
            m = max(m, float(np.max(np.abs(C.periodic_xcorr(codes[i], codes[j])))))
    return m


CROSS = {'m-sequence, N=127': cross_max(mseq7),
         'Gold, N=127': cross_max(gold7),
         'Kasami small, N=63': cross_max(kas6),
         'Golay pair, N=128 (2 shots)': float(np.max(np.abs(C.periodic_xcorr(ga, gb))))}


def run(K, nchips, ac, rng):
    """RMS Bragg error with code leakage and detector noise only."""
    bins = np.sort(rng.choice(np.arange(1, nchips), size=K, replace=False))
    nub = rng.uniform(-DETUNE_GHZ, DETUNE_GHZ, size=K)
    A = R * np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
    W = ac[(bins[:, None] - bins[None, :]) % nchips]
    np.fill_diagonal(W, 0.0)
    S = A + W @ A + rng.normal(0, SIGMA_N, size=A.shape)
    errs = np.array([(C.gauss_fit_peak(nu, S[k]) - nub[k]) * PM for k in range(K)])
    return float(np.sqrt(np.mean(errs ** 2)))


Ks = np.array([4, 8, 16, 24, 32, 48])
NT = 8
curves = {}
for name, nchips, fsize, ac in FAMILIES:
    vals = []
    for K in Ks:
        if K >= nchips:
            vals.append(np.nan); continue
        vals.append(np.mean([run(K, nchips, ac, np.random.default_rng(5000 + K + t))
                             for t in range(NT)]))
    curves[name] = np.array(vals)

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
COL = {'m-sequence, N=127': '#2980b9', 'Gold, N=127': '#c0392b',
       'Kasami small, N=63': '#7d3c98', 'Golay pair, N=128 (2 shots)': '#f39c12'}
fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.7))
for name, nchips, fsize, ac in FAMILIES:
    L = len(ac); lag = np.arange(L) - L // 2
    y = np.abs(np.roll(ac, L // 2))
    sel = np.abs(lag) <= 30
    ax[0].semilogy(lag[sel], np.maximum(y[sel], 1e-4), lw=1.1, color=COL[name], label=name)
ax[0].axhline(17 / 127, color='0.5', ls=':', lw=0.9)
ax[0].text(30, 17 / 127 * 1.3, 'Gold bound 17/127 ', fontsize=7, color='0.4', ha='right')
ax[0].axhline(1 / 127, color='0.5', ls='--', lw=0.9)
ax[0].text(30, 1 / 127 * 1.3, 'm-sequence 1/N ', fontsize=7, color='0.4', ha='right')
ax[0].set_ylim(1e-4, 2)
ax[0].set_xlabel('lag [chip]'); ax[0].set_ylabel('|periodic autocorrelation|')
ax[0].set_title('(a) Autocorrelation: what separates gratings')
ax[0].legend(fontsize=7, loc='lower center'); ax[0].grid(True, which='both', alpha=0.25)

for name, nchips, fsize, ac in FAMILIES:
    ax[1].plot(Ks, curves[name], 'o-', color=COL[name], label=name)
ax[1].axhline(10.0, color='0.3', ls='--', lw=0.8)
ax[1].text(0.97, 0.06, '10 pm target', transform=ax[1].transAxes, fontsize=7.5,
           color='0.3', ha='right')
ax[1].set_xlabel('gratings on the fiber, K')
ax[1].set_ylabel('RMS Bragg error [pm]')
ax[1].set_title('(b) Code-induced error vs array size')
ax[1].legend(fontsize=7, loc='upper left'); ax[1].grid(True, alpha=0.25)
plt.tight_layout(); plt.savefig('figs/fig_s14_codes.png', dpi=140)

# ---------------------------------------------------------------------------
# printed tables
# ---------------------------------------------------------------------------
print('code families for a swept CDM interrogator (R=%.2f, M=%d slots)' % (R, M_SLOTS))
print('--- properties ---')
print('   family                          N   size   max |auto sidelobe|   max |cross|')
for name, nchips, fsize, ac in FAMILIES:
    side = float(np.max(np.abs(ac[1:])))
    print('   %-28s %4d   %4d   %19.4f   %11.4f'
          % (name, nchips, fsize, side, CROSS[name]))
print('--- (b) RMS Bragg error [pm] from code leakage alone ---')
print('   family                        ' + ''.join('  K=%-5d' % K for K in Ks))
for name, nchips, fsize, ac in FAMILIES:
    row = ''.join('  %6.2f ' % v if np.isfinite(v) else '     -  ' for v in curves[name])
    print('   %-28s%s' % (name, row))
i32 = int(np.where(Ks == 32)[0][0])
print('  => at K = 32 the m-sequence beats Gold by %.1fx (%.2f vs %.2f pm)'
      % (curves['Gold, N=127'][i32] / curves['m-sequence, N=127'][i32],
         curves['m-sequence, N=127'][i32], curves['Gold, N=127'][i32]))
print('  => the Golay pair removes the side lobe entirely, at twice the acquisition time')
print('  => Gold only pays off if several wavelengths are simultaneously on the fiber,'
      ' which a swept source avoids')
print('saved figs/fig_s14_codes.png')
