"""s22_codelength.py - the code length as a design variable.

The sequence length N is the one number the designer types into the signal
generator, and it buys four different things at once. This script puts all four
on one chart so the trade can be made deliberately.

  1. DELAY BINS AND RANGE. N chips give N resolvable delay bins and an
     unambiguous range of N * c/(2nB); beyond it the periodic correlation wraps.

  2. LEAKAGE CEILING. The side-lobe floor of an m-sequence is 1/N, so the
     leakage-limited array size grows linearly with N (the K = 0.37 N rule of the
     capacity study). Verified here for four lengths.

  3. MEASUREMENT TIME. One sequence period is N/B; equivalent-time sampling
     needs a few periods per wavelength step, so the frame time and refresh rate
     fall as N grows. Long codes are not free.

  4. COLLISION STATISTICS. With K gratings dropped on N bins, the number of
     co-binned pairs follows the birthday problem: the expected count is
     K(K-1)/(2N). A co-binned pair is not lost: the two-component fit of the
     inversion study resolves it whenever the two Bragg wavelengths differ by at
     least half a linewidth. Only the remaining fraction, co-binned AND closer
     than half a linewidth in wavelength, needs anything beyond a parametric
     fit, and that is the honest size of the niche for learned demodulators in
     this architecture. For detunings uniform over the +-200 pm Peltier window
     that fraction is P_hard = 1 - (1 - d/W)^2 with d = FWHM/2 = 125 pm and
     W = 400 pm, about 0.53, so the expected number of hard pairs is roughly
     0.26 K^2 / N. For K = 32 on N = 511 that is half a pair per array.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS
FS.apply()


PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
C_LIGHT, N_GROUP = 2.99792458e8, 1.468
M_SLOTS = 64
DETUNE = 25.0
nu = np.linspace(-2.6 * F, 2.6 * F, M_SLOTS)
B = 25e6
FS = 20e6
OSR = 4
R = 0.10
SIGMA_N = 1.1e-6

NBITS = [6, 7, 8, 9, 10]
NS = [2 ** n - 1 for n in NBITS]
# number of m-sequences of degree n: euler_phi(2^n - 1) / n
FAMILY = {63: 6, 127: 18, 255: 16, 511: 48, 1023: 60}


def acorr_of(nbits):
    m = 1.0 - 2.0 * C._mls01(nbits)
    return C.periodic_xcorr(m, m)


def leak_error(K, nchips, ac, seed=0, ntrials=4):
    """Leakage plus noise only, as in the capacity study."""
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


def k_at_target(nchips, ac, target=10.0):
    Ks = np.array([4, 8, 16, 24, 32, 48, 64, 96, 128, 192])
    Ks = Ks[Ks < nchips]
    e = np.array([leak_error(K, nchips, ac, seed=1000 + K) for K in Ks])
    if e[-1] <= target:
        return float(Ks[-1]), Ks, e, True
    i = int(np.argmax(e > target))
    if i == 0:
        return 0.0, Ks, e, False
    x0, x1, y0, y1 = Ks[i - 1], Ks[i], e[i - 1], e[i]
    return float(x0 + (target - y0) * (x1 - x0) / (y1 - y0)), Ks, e, False


def ets_step(nchips, chip_rate=B, fs=FS, osr=OSR):
    t_seq = nchips / chip_rate
    per = max(1, int(np.floor(fs * t_seq)))
    periods = int(np.ceil(nchips * osr / per))
    return periods * t_seq


kmax, sat, curves = {}, {}, {}
for n, N in zip(NBITS, NS):
    ac = acorr_of(n)
    kmax[N], Ks_, e_, sat[N] = k_at_target(N, ac)
    curves[N] = (Ks_, e_)

# collision statistics
P_HARD = 1.0 - (1.0 - 125.0 / 400.0) ** 2      # |dLambda| < FWHM/2, triangular


def p_collision(K, N):
    i = np.arange(1, K)
    return 1.0 - np.prod(1.0 - i / N)


Kx = np.arange(2, 130)

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
# single panel: the co-binning statistics moved to eq. (18) in the text
fig, ax0 = plt.subplots(figsize=(3.45, 2.6))
ax = [ax0]

NSarr = np.array(NS, float)
km = np.array([kmax[N] for N in NS])
solid = np.array([not sat[N] for N in NS])
ax[0].loglog(NSarr[solid], km[solid], 'o-', color='#2980b9',
             label='simulated, 10 pm target')
ax[0].loglog(NSarr[~solid], km[~solid], 'o', mfc='none', mec='#2980b9',
             label='lower bound (sweep cap)')
ax[0].loglog(NSarr, 0.37 * NSarr, '--', color='0.4', label=r'$K = 0.37\,N$')
for N in NS:
    lab = 'N=%d' % N
    dx, dy = (4, -11) if N < 1023 else (-30, 8)
    ax[0].annotate(lab, (N, kmax[N]), textcoords='offset points',
                   xytext=(dx, dy), fontsize=6.6, color='#2980b9')
axr = ax[0].twinx()
rng_km = NSarr * C_LIGHT / (2 * N_GROUP * B) / 1000.0
axr.loglog(NSarr, rng_km, 's:', color='#c0392b',
           label='alias-free reach $Nc/(2nB)$')
axr.set_ylabel('alias-free reach $Nc/(2nB)$ [km]\n(code period caps the'
               ' delay window)', color='#c0392b', fontsize=7)
axr.tick_params(axis='y', labelcolor='#c0392b', labelsize=7)
ax[0].set_xlabel('code length N')
ax[0].set_ylabel('leakage-limited gratings per band', color='#2980b9')
ax[0].tick_params(axis='y', labelcolor='#2980b9', labelsize=7)
ax[0].set_xticks(NS); ax[0].set_xticklabels([str(N) for N in NS], fontsize=7)
ax[0].legend(fontsize=6.4, loc='upper left')
ax[0].grid(True, which='both', alpha=0.25)

fig.tight_layout()
fig.savefig('figs/fig_s22_codelength.png', dpi=150, bbox_inches='tight')
fig.savefig('figs/fig_s22_codelength.pdf', bbox_inches='tight')

# ---------------------------------------------------------------------------
# design table
# ---------------------------------------------------------------------------
print('CODE LENGTH AS A DESIGN VARIABLE  (B = %.0f Mchip/s, ADC %.0f MS/s,'
      ' OSR = %d, M = %d wavelength steps)' % (B / 1e6, FS / 1e6, OSR, M_SLOTS))
print('     N   family   T_seq[us]   range[m]   K_leak(10pm)   step[us]'
      '   frame[ms]   refresh[Hz]')
for N in NS:
    t_seq = N / B * 1e6
    rng_m = N * C_LIGHT / (2 * N_GROUP * B)
    step = ets_step(N) * 1e6
    frame = M_SLOTS * step / 1e3
    ktxt = ('>= %.0f' % kmax[N]) if sat[N] else ('%.1f' % kmax[N])
    print('  %5d   %6d   %9.2f   %8.0f   %12s   %8.1f   %9.2f   %11.0f'
          % (N, FAMILY[N], t_seq, rng_m, ktxt, step, frame, 1000.0 / frame))
print()
print('collision statistics (random placement):')
print('  expected co-binned pairs = K(K-1)/2N; hard fraction (closer than'
      ' FWHM/2) = %.2f' % P_HARD)
for K, N in [(32, 127), (32, 511), (64, 511), (64, 1023)]:
    pairs = K * (K - 1) / (2.0 * N)
    print('  K=%3d, N=%4d: P(any shared bin) = %.2f, expected pairs = %.2f,'
          ' expected hard pairs = %.2f'
          % (K, N, p_collision(K, N), pairs, P_HARD * pairs))
print('saved figs/fig_s22_codelength.png')
