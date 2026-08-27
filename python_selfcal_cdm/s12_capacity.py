"""s12_capacity.py - how many gratings fit on one code-multiplexed fiber?

Architecture (fixed by the project lead): the spreading CODE marks the optical
wavelength, the round-trip DELAY marks the grating position. The source is a
swept VCSEL, so wavelength slots are visited SEQUENTIALLY: inside one correlation
window a single code is on the fiber and every grating returns that same code,
displaced by its own round-trip delay. Two consequences set the error model:

  * gratings are separated by the code AUTOCORRELATION side lobes, not by
    cross-correlation between different codes. A maximal-length sequence has a
    side lobe of exactly -1/N; a Gold code is three valued and reaches
    (2^((n+1)/2)+1)/N, i.e. 17/127 here. Cross-correlation only matters when
    several wavelengths are simultaneously on the fiber, which a swept source
    avoids. This reverses the usual "Gold is better" reflex;
  * the leakage in one delay bin is a weighted sum of the other gratings'
    spectra, so it is a broad, smooth background under a narrow peak, and a fit
    with a baseline term removes most of it.

Mechanisms, all evaluated in the despread domain (despreading is linear, so the
correlator is replaced by its output, while the real measured autocorrelation of
the actual code supplies the leakage weights):

  1. SPECTRAL SHADOWING. Grating k is illuminated through the upstream gratings,
     giving prod_{j<k} (1 - R*S_j(nu))^2. With gratings that share a nominal
     wavelength this notch is deep and, when the upstream gratings are detuned,
     asymmetric across the line of grating k, so it biases the fitted peak. The
     effect is cumulative, which is what limits strong-grating arrays.

  2. MULTIPLE REFLECTIONS (ghosts). A third-order path reflecting at gratings
     a, b, c (b upstream of both a and c) returns at apparent delay
     tau_a - tau_b + tau_c with power ~ R^3. The array is incoherent over its
     length (s9: L_coh = 0.22 m << spacing), so powers add and the ghost line is
     the PRODUCT of three grating power spectra: narrower than a real peak, as
     reported by Guo 2017, and centred at the mean of the three detunings.

     Structural result: for UNIFORM spacing tau_k = tau_0 + k*D the ghost delay
     is tau_0 + D*(a-b+c), which is again a multiple of D, so every ghost that
     falls inside the span of the array lands exactly on an occupied bin.
     Randomised spacing scatters almost all of them into empty bins, where they
     are harmless.

  3. CODE LEAKAGE and DETECTOR NOISE. Noise comes from an explicit link budget
     (source power, circulator and connector loss, APD noise-equivalent power,
     chip-rate bandwidth, correlation processing gain, unipolar modulation),
     not from an assumed SNR, so weak gratings are penalised only as much as the
     physics penalises them.

Outputs
  (a) RMS Bragg error vs array size at R = 10 percent, decomposed by mechanism.
  (b) Capacity (largest array meeting a 10 pm target) vs reflectivity, for
      uniform and randomised spacing.
"""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS
FS.apply()

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ            # 31.25 GHz = 250 pm (procured FBGS gratings)
SIG = F / 2.35482

# --- design parameters held fixed in this study -----------------------------
NBITS = 7
N_CHIPS = 2 ** NBITS - 1      # 127 chips -> 127 resolvable delay bins
M_SLOTS = 64                  # wavelength samples across the sweep
BAND_GHZ = 2.6 * F            # sweep half-span
DETUNE_GHZ = 25.0             # spread of grating detunings (temperature/strain)
TARGET_PM = 10.0              # accuracy target defining "capacity"

# --- optical link budget (sets the absolute noise floor) --------------------
P_SOURCE_W = 1.0e-3           # BW10-1550 fiber-coupled output, order of 0 dBm
LOSS_DB = 4.0                 # circulator round trip plus connectors
NEP_W_RTHZ = 0.5e-12          # InGaAs APD module noise-equivalent power
CHIP_RATE_HZ = 25.0e6         # sets the detection bandwidth
UNIPOLAR = 0.5                # on/off keying recovers half the bipolar amplitude
_p_rx = P_SOURCE_W * 10.0 ** (-LOSS_DB / 10.0)
_p_noise = NEP_W_RTHZ * np.sqrt(CHIP_RATE_HZ)
SIGMA_N = (_p_noise / _p_rx) / (np.sqrt(N_CHIPS) * UNIPOLAR)   # referred to R = 1

nu = np.linspace(-BAND_GHZ, BAND_GHZ, M_SLOTS)

# real code autocorrelation, used as the leakage weight between delay bins
_MSEQ = 1.0 - 2.0 * C._mls01(NBITS)
ACORR_MSEQ = C.periodic_xcorr(_MSEQ, _MSEQ)
_GOLD = C.gold_set(NBITS)[3]
ACORR_GOLD = C.periodic_xcorr(_GOLD, _GOLD)


def _layout(K, mode, rng, nchips=None):
    """Delay bins (integer chips) for K gratings inside nchips bins."""
    nchips = N_CHIPS if nchips is None else nchips
    if mode == 'uniform':
        step = max(1, (nchips - 2) // K)
        return 1 + step * np.arange(K)
    return np.sort(rng.choice(np.arange(1, nchips), size=K, replace=False))


def _spectra(K, R, rng, with_shadow=True):
    nub = rng.uniform(-DETUNE_GHZ, DETUNE_GHZ, size=K)
    shapes = np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)   # [K, M]
    tcum = np.ones((K, M_SLOTS))
    if with_shadow:
        trans = (1.0 - R * shapes) ** 2
        for k in range(1, K):
            tcum[k] = tcum[k - 1] * trans[k - 1]
    return nub, shapes, tcum


def _ghost_matrix(K, nub, R, bins, tcum, chunk=40000):
    """Sum of third-order ghosts landing in each occupied delay bin -> [K, M]."""
    idx = np.arange(K)
    a, b, c = np.meshgrid(idx, idx, idx, indexing='ij')
    m = (b < a) & (b < c)
    a, b, c = a[m], b[m], c[m]
    gbin = bins[a] - bins[b] + bins[c]
    order = np.argsort(bins)
    pos = np.searchsorted(bins[order], gbin)
    pos = np.clip(pos, 0, K - 1)
    hit = bins[order][pos] == gbin                    # ghost lands on a real bin
    if not hit.any():
        return np.zeros((K, M_SLOTS))
    a, b, c = a[hit], b[hit], c[hit]
    target = order[pos[hit]]                          # which grating it pollutes
    cbar = (nub[a] + nub[b] + nub[c]) / 3.0
    spread = (nub[a] - cbar) ** 2 + (nub[b] - cbar) ** 2 + (nub[c] - cbar) ** 2
    amp = R ** 3 * np.exp(-0.5 * spread / SIG ** 2)   # product of three Gaussians
    upstream = np.minimum(np.minimum(a, b), c)
    sig_g = SIG / np.sqrt(3.0)
    out = np.zeros((K, M_SLOTS))
    for i0 in range(0, len(amp), chunk):
        s = slice(i0, i0 + chunk)
        g = (amp[s, None] * tcum[upstream[s]] *
             np.exp(-0.5 * ((nu[None, :] - cbar[s, None]) / sig_g) ** 2))
        np.add.at(out, target[s], g)
    return out


def run(K, R, mode, rng, acorr=None, with_ghosts=True, with_shadow=True,
        with_leak=True, with_noise=True, nchips=None, peel=False):
    """RMS Bragg-position error over a K-grating array [pm]."""
    nchips = N_CHIPS if nchips is None else nchips
    acorr = ACORR_MSEQ if acorr is None else acorr
    bins = _layout(K, mode, rng, nchips)
    nub, shapes, tcum = _spectra(K, R, rng, with_shadow=with_shadow)
    primary = R * shapes * tcum                                     # [K, M]
    G = _ghost_matrix(K, nub, R, bins, tcum) if with_ghosts else 0.0

    if with_leak:
        W = acorr[(bins[:, None] - bins[None, :]) % nchips]
        np.fill_diagonal(W, 0.0)
        leak = W @ primary
    else:
        leak = 0.0
    S = primary + leak + G
    if with_noise:
        S = S + rng.normal(0, SIGMA_N, size=S.shape)
    if peel:
        errs = _peel(S, nub)
    else:
        errs = np.array([(C.gauss_fit_peak(nu, S[k]) - nub[k]) * PM for k in range(K)])
    return float(np.sqrt(np.mean(errs ** 2)))


def _peel(S, nub, floor=0.05):
    """Sequential deshadowing.

    Shadowing is deterministic once the upstream gratings are known, so it can be
    undone rather than merely tolerated. Grating 1 is unshadowed, so its line is
    fitted directly; from the fitted amplitude, centre and width the transmission
    it imposes is reconstructed and divided out of every grating behind it. The
    procedure then repeats down the array. Only the array ORDER is needed, which
    the delay bins already give."""
    K = S.shape[0]
    tcorr = np.ones(M_SLOTS)
    errs = np.empty(K)
    for k in range(K):
        Sk = S[k] / np.maximum(tcorr, floor)
        a, mu, sg, _ = C.gauss_fit_full(nu, Sk)
        errs[k] = (mu - nub[k]) * PM
        line = np.clip(a, 0.0, 0.99) * np.exp(-0.5 * ((nu - mu) / max(sg, 1e-3)) ** 2)
        tcorr = tcorr * (1.0 - line) ** 2
    return errs


def sweep(Ks, R, mode, ntrials, seed, **kw):
    out = []
    for K in Ks:
        rng = np.random.default_rng(seed + K)
        out.append(np.mean([run(K, R, mode, rng, **kw) for _ in range(ntrials)]))
    return np.array(out)


# ---------------------------------------------------------------------------
# (a) error vs array size at R = 10 percent, decomposed
# ---------------------------------------------------------------------------
Ks = np.array([4, 8, 16, 24, 32, 48, 64, 96])
NT = 4
R_A = 0.10
full_rnd = sweep(Ks, R_A, 'random', NT, 100)
full_uni = sweep(Ks, R_A, 'uniform', NT, 100)
only_shadow = sweep(Ks, R_A, 'random', NT, 100, with_ghosts=False, with_leak=False,
                    with_noise=False)
only_ghost = sweep(Ks, R_A, 'random', NT, 100, with_shadow=False, with_leak=False,
                   with_noise=False)
only_ghost_uni = sweep(Ks, R_A, 'uniform', NT, 100, with_shadow=False, with_leak=False,
                       with_noise=False)
only_leak = sweep(Ks, R_A, 'random', NT, 100, with_shadow=False, with_ghosts=False)
peeled_rnd = sweep(Ks, R_A, 'random', NT, 100, peel=True)

# ---------------------------------------------------------------------------
# (b) capacity vs reflectivity
# ---------------------------------------------------------------------------
Rs = np.array([0.0003, 0.001, 0.003, 0.01, 0.03, 0.05, 0.10, 0.20, 0.30])
Ks_cap = np.array([4, 8, 16, 32, 48, 64, 96])


def capacity(R, mode, seed, **kw):
    e = sweep(Ks_cap, R, mode, 3, seed, **kw)
    if e[0] > TARGET_PM:
        return 0.0
    if e[-1] <= TARGET_PM:
        return float(Ks_cap[-1])
    i = int(np.argmax(e > TARGET_PM))
    x0, x1 = Ks_cap[i - 1], Ks_cap[i]
    y0, y1 = e[i - 1], e[i]
    return float(x0 + (TARGET_PM - y0) * (x1 - x0) / (y1 - y0))


cap_rnd = np.array([capacity(R, 'random', 500 + i * 37) for i, R in enumerate(Rs)])
cap_uni = np.array([capacity(R, 'uniform', 900 + i * 37) for i, R in enumerate(Rs)])
cap_peel = np.array([capacity(R, 'random', 500 + i * 37, peel=True) for i, R in enumerate(Rs)])

# code family check at a common operating point
K_CHK, R_CHK = 32, 0.03
e_mseq = np.mean([run(K_CHK, R_CHK, 'random', np.random.default_rng(11 + t),
                      acorr=ACORR_MSEQ) for t in range(5)])
e_gold = np.mean([run(K_CHK, R_CHK, 'random', np.random.default_rng(11 + t),
                      acorr=ACORR_GOLD) for t in range(5)])

# bench configuration: 3 procured gratings, R = 10 percent
bench = np.mean([run(3, 0.10, 'random', np.random.default_rng(700 + t)) for t in range(8)])
bench_ns = np.mean([run(3, 0.10, 'random', np.random.default_rng(700 + t),
                        with_shadow=False) for t in range(8)])
bench_pl = np.mean([run(3, 0.10, 'random', np.random.default_rng(700 + t),
                        peel=True) for t in range(8)])

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
FLOOR = 0.2   # plotting floor: curves below this are not resolvable anyway
cl = lambda a: np.maximum(a, FLOOR)
fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.7))
ax[0].semilogy(Ks, cl(only_shadow), '^--', color='#7d3c98', label='spectral shadowing only')
ax[0].semilogy(Ks, cl(only_ghost_uni), 'v--', color='#c0392b', label='ghosts only, uniform spacing')
ax[0].semilogy(Ks, cl(only_ghost), 'v:', color='#28b463', label='ghosts only, randomised spacing')
ax[0].semilogy(Ks, cl(only_leak), 'd--', color='0.55', label='code leakage + noise only')
ax[0].semilogy(Ks, cl(full_rnd), 's-', color='#2980b9', lw=1.8, label='full model, randomised')
ax[0].semilogy(Ks, cl(peeled_rnd), 'o-', color='#f39c12', lw=1.8,
               label='full model + sequential deshadowing')
ax[0].axhline(TARGET_PM, color='0.3', ls=':', lw=1.0)
ax[0].text(Ks[-1], TARGET_PM * 1.25, '%.0f pm target ' % TARGET_PM, fontsize=7.5,
           color='0.3', ha='right')
ax[0].set_ylim(FLOOR, 400)
ax[0].set_xlabel('gratings on the fiber, K')
ax[0].set_ylabel('RMS Bragg error [pm]')
ax[0].set_title('(a) Error mechanisms vs array size, R = %.0f%%' % (R_A * 100))
ax[0].legend(fontsize=6.8, loc='lower right'); ax[0].grid(True, which='both', alpha=0.25)

ax[1].semilogx(Rs * 100, cap_uni, 'o-', color='#c0392b', label='uniform spacing')
ax[1].semilogx(Rs * 100, cap_rnd, 's-', color='#2980b9', label='randomised spacing')
ax[1].semilogx(Rs * 100, cap_peel, '^-', color='#f39c12', label='randomised + deshadowing')
ax[1].axvline(0.10 * 100, color='0.5', ls=':', lw=1.0)
ax[1].text(0.10 * 100, 52, 'procured R = 10% ', fontsize=7, ha='right', color='0.3')
ax[1].set_ylim(0, 58)
ax[1].set_xlabel('grating reflectivity R [%]')
ax[1].set_ylabel('gratings meeting the %.0f pm target' % TARGET_PM)
ax[1].set_title('(b) Capacity vs reflectivity')
ax[1].legend(fontsize=7.5, loc='lower left'); ax[1].grid(True, which='both', alpha=0.25)
plt.tight_layout(); plt.savefig('figs/fig_s12_capacity.png', dpi=140)
plt.savefig('figs/fig_s12_capacity.pdf')

# ---------------------------------------------------------------------------
# printed tables
# ---------------------------------------------------------------------------
print('CDM capacity study: N=%d chips, M=%d wavelength slots, FWHM=%.0f pm,'
      ' detuning spread +-%.0f GHz' % (N_CHIPS, M_SLOTS, C.FBG_FWHM_PM, DETUNE_GHZ))
print('link budget: %.0f dBm source, %.0f dB loss, NEP %.1f pW/rtHz, B = %.0f Mchip/s,'
      % (10 * np.log10(P_SOURCE_W / 1e-3), LOSS_DB, NEP_W_RTHZ * 1e12, CHIP_RATE_HZ / 1e6))
print('             processing gain sqrt(%d), unipolar factor %.1f -> sigma = %.2e (referred to R=1)'
      % (N_CHIPS, UNIPOLAR, SIGMA_N))
print('             equivalent SNR of a single R=10%% grating: %.0f dB'
      % (20 * np.log10(0.10 / SIGMA_N)))
print('--- (a) RMS Bragg error [pm] vs array size, R = %.0f%% ---' % (R_A * 100))
print('     K   full/rand   shadowing   ghost/uni  ghost/rand  leak+noise   deshadowed')
for i, K in enumerate(Ks):
    print('  %4d  %10.2f  %10.2f  %10.2f  %10.2f  %10.2f  %11.2f'
          % (K, full_rnd[i], only_shadow[i], only_ghost_uni[i],
             only_ghost[i], only_leak[i], peeled_rnd[i]))
print('--- (b) capacity at the %.0f pm target ---' % TARGET_PM)
print('     R[%]   K_max random   K_max uniform   K_max deshadowed')
for i, R in enumerate(Rs):
    print('  %7.2f   %12.1f   %13.1f   %16.1f'
          % (R * 100, cap_rnd[i], cap_uni[i], cap_peel[i]))
ibest = int(np.argmax(cap_rnd))
print('  => randomised spacing peaks at R = %.2f%% with K_max = %.0f'
      % (Rs[ibest] * 100, cap_rnd[ibest]))
print('--- code family, K=%d, R=%.0f%%, randomised spacing ---' % (K_CHK, R_CHK * 100))
print('  m-sequence (side lobe 1/N)  : %.2f pm' % e_mseq)
print('  Gold code (three valued)    : %.2f pm' % e_gold)
# leakage ceiling vs code length: leakage/primary is independent of R, so this
# ceiling holds whatever the reflectivity
_M9 = 1.0 - 2.0 * C._mls01(9)
ACORR_M9 = C.periodic_xcorr(_M9, _M9)
print('--- leakage-only ceiling vs code length (R = 0.3%%, randomised spacing) ---')
print('       N      K    RMS error [pm]')
for nch, ac, tag in [(127, ACORR_MSEQ, '127'), (511, ACORR_M9, '511')]:
    for K in [16, 32, 48, 64, 96]:
        e = np.mean([run(K, 0.003, 'random', np.random.default_rng(3100 + t + K),
                         acorr=ac, nchips=nch, with_shadow=False, with_ghosts=False)
                     for t in range(3)])
        print('  %6s   %4d   %14.2f' % (tag, K, e))
print('--- bench configuration: 3 procured gratings, R = 10%% ---')
print('  full model            : %.2f pm' % bench)
print('  shadowing switched off: %.2f pm' % bench_ns)
print('  sequential deshadowing: %.2f pm' % bench_pl)
print('saved figs/fig_s12_capacity.png')
