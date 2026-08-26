"""s13_chiprate.py - what the acquisition chain must deliver.

The correlator separates gratings by round-trip delay, so the chip rate B fixes
the spatial resolution dz = c / (2 n B): 25 Mchip/s gives 4.1 m in fiber,
50 Mchip/s gives 2.0 m. Gratings closer than dz share a delay bin, their spectra
add, and the fitted Bragg position of each is pulled toward the other. This
script turns that statement into numbers and then asks what the electronics must
provide: timing resolution, timing jitter, converter resolution, and how long an
equivalent-time acquisition takes.

Model. The despread readout of the bin belonging to grating k is the sum of every
grating's spectrum weighted by the code autocorrelation evaluated at the delay
difference. Inside one chip the autocorrelation of a rectangular-chip sequence is
the triangle 1 - |dt| / T_chip; beyond it, the periodic side lobe -1/N. Timing
jitter is applied as a Gaussian blur of that triangle, which both lowers the peak
and widens the skirt, so a neighbour leaks more.

Outputs
  (a) RMS Bragg error vs chip rate for three minimum grating spacings, with the
      resolution criterion dz = d_min marked.
  (b) RMS Bragg error vs RMS timing jitter, in units of the chip period.
  Printed: equivalent-time sampling budget and converter-resolution budget.
"""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import common as C

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
C_LIGHT = 2.99792458e8
N_GROUP = 1.468                      # SMF-28 group index

NBITS = 7
N_CHIPS = 2 ** NBITS - 1
M_SLOTS = 64
BAND_GHZ = 2.6 * F
DETUNE_GHZ = 25.0
R = 0.10
K = 8                                # gratings in the test array
SIDELOBE = -1.0 / N_CHIPS            # m-sequence periodic side lobe

nu = np.linspace(-BAND_GHZ, BAND_GHZ, M_SLOTS)


def dz_of(B):
    """Spatial resolution in fiber for chip rate B [Hz], in metres."""
    return C_LIGHT / (2.0 * N_GROUP * B)


def weights(dt_chips, jitter_chips=0.0):
    """Code autocorrelation sampled at a delay difference in chips, with the
    triangular main lobe optionally blurred by RMS timing jitter."""
    if jitter_chips <= 0:
        main = np.clip(1.0 - np.abs(dt_chips), 0.0, None)
    else:
        u = np.linspace(-4 * jitter_chips, 4 * jitter_chips, 41)
        g = np.exp(-0.5 * (u / jitter_chips) ** 2); g /= g.sum()
        main = np.zeros_like(dt_chips, dtype=float)
        for ui, gi in zip(u, g):
            main += gi * np.clip(1.0 - np.abs(dt_chips - ui), 0.0, None)
    return np.where(np.abs(dt_chips) < 1.0 + 4 * jitter_chips, main, 0.0) + \
        SIDELOBE * (np.abs(dt_chips) >= 1.0 + 4 * jitter_chips)


def layout(rng, d_min, span=40.0):
    """Grating positions [m] with a guaranteed minimum spacing."""
    z = [2.0]
    while len(z) < K:
        cand = z[-1] + d_min * rng.uniform(1.0, 1.0 + span / (K * d_min))
        z.append(cand)
    return np.array(z)


def run(B, d_min, rng, jitter_chips=0.0):
    z = layout(rng, d_min)
    tau = 2.0 * N_GROUP * z / C_LIGHT            # round-trip delay [s]
    t_chips = tau * B
    nub = rng.uniform(-DETUNE_GHZ, DETUNE_GHZ, size=K)
    A = R * np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
    errs = np.empty(K)
    for k in range(K):
        lag = np.round(t_chips[k])
        w = weights(lag - t_chips, jitter_chips=jitter_chips)
        S = w @ A
        errs[k] = (C.gauss_fit_peak(nu, S) - nub[k]) * PM
    return float(np.sqrt(np.mean(errs ** 2)))


def sweep(Bs, d_min, ntrials=12, seed=0, jitter_chips=0.0):
    out = []
    for B in Bs:
        rng = np.random.default_rng(seed + int(B / 1e5))
        out.append(np.mean([run(B, d_min, rng, jitter_chips) for _ in range(ntrials)]))
    return np.array(out)


# ---------------------------------------------------------------------------
# (a) error vs chip rate for three minimum spacings
# ---------------------------------------------------------------------------
Bs = np.array([2, 3, 5, 8, 12, 18, 25, 35, 50, 75, 100, 150, 200]) * 1e6
e_05 = sweep(Bs, 0.5, seed=10)
e_2 = sweep(Bs, 2.0, seed=20)
e_4 = sweep(Bs, 4.0, seed=30)

# ---------------------------------------------------------------------------
# (b) error vs timing jitter at a working point
# ---------------------------------------------------------------------------
B_WORK = 25e6
jit = np.array([0.0, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50])
e_jit_25 = np.array([sweep(np.array([25e6]), 4.0, ntrials=14, seed=40,
                           jitter_chips=j)[0] for j in jit])
e_jit_50 = np.array([sweep(np.array([50e6]), 4.0, ntrials=14, seed=40,
                           jitter_chips=j)[0] for j in jit])

# ---------------------------------------------------------------------------
# equivalent-time sampling budget
# ---------------------------------------------------------------------------
OSR = 4                              # samples per chip in the reconstructed record


def ets(B, fs):
    """Equivalent-time acquisition of one full record.

    The probe is periodic with T_seq = N / B. An ADC running at fs takes
    floor(fs*T_seq) samples per period; the sampling instant is stepped by
    1/(OSR*B) between periods until N*OSR points are collected."""
    t_seq = N_CHIPS / B
    per_period = max(1, int(np.floor(fs * t_seq)))
    n_total = N_CHIPS * OSR
    periods = int(np.ceil(n_total / per_period))
    return periods * t_seq, per_period, periods, 1.0 / (OSR * B)


# ---------------------------------------------------------------------------
# converter resolution budget
# ---------------------------------------------------------------------------
def quant_jitter(bits, headroom=0.5, ntrials=300, seed=7):
    """RMS Bragg jitter from ADC quantisation alone.

    The grating peak uses `headroom` of full scale; quantisation noise is
    FS / (2^bits * sqrt(12)) and the correlator gives sqrt(N) of processing gain."""
    rng = np.random.default_rng(seed)
    fs_units = R / headroom
    sigma = fs_units / (2.0 ** bits * np.sqrt(12.0)) / np.sqrt(N_CHIPS)
    errs = []
    for _ in range(ntrials):
        nub = rng.uniform(-DETUNE_GHZ, DETUNE_GHZ)
        S = R * np.exp(-0.5 * ((nu - nub) / SIG) ** 2) + rng.normal(0, sigma, M_SLOTS)
        errs.append((C.gauss_fit_peak(nu, S) - nub) * PM)
    return float(np.std(errs)), sigma


# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
plt.rcParams.update({'font.size': 9})
fig, ax = plt.subplots(1, 2, figsize=(10.6, 4.0))
ax[0].loglog(Bs / 1e6, np.maximum(e_05, 0.05), 'o-', color='#c0392b',
             label='min spacing 0.5 m')
ax[0].loglog(Bs / 1e6, np.maximum(e_2, 0.05), 's-', color='#2980b9',
             label='min spacing 2 m')
ax[0].loglog(Bs / 1e6, np.maximum(e_4, 0.05), '^-', color='#28b463',
             label='min spacing 4 m (bench)')
for d, col in [(0.5, '#c0392b'), (2.0, '#2980b9'), (4.0, '#28b463')]:
    B_req = C_LIGHT / (2.0 * N_GROUP * d) / 1e6
    ax[0].axvline(B_req, color=col, ls=':', lw=0.9)
ax[0].axhline(10.0, color='0.3', ls='--', lw=0.8)
ax[0].text(0.03, 0.10, 'dotted lines: dz = min spacing', transform=ax[0].transAxes,
           fontsize=7, color='0.3')
ax[0].text(0.03, 0.04, 'dashed line: 10 pm target', transform=ax[0].transAxes,
           fontsize=7, color='0.3')
ax[0].set_xlabel('chip rate B [Mchip/s]')
ax[0].set_ylabel('RMS Bragg error [pm]')
ax[0].set_title('(a) Resolution: error vs chip rate')
ax[0].legend(fontsize=7.5, loc='upper right'); ax[0].grid(True, which='both', alpha=0.25)

ax[1].plot(jit, np.maximum(e_jit_25, 0.05), '^-', color='#28b463',
               label='B = 25 Mchip/s (dz = 4.1 m, marginal)')
ax[1].plot(jit, np.maximum(e_jit_50, 0.05), 's-', color='#2980b9',
               label='B = 50 Mchip/s (dz = 2.0 m, comfortable)')
ax[1].axhline(10.0, color='0.3', ls='--', lw=0.8)
ax[1].text(0.97, 0.93, '10 pm target', transform=ax[1].transAxes, fontsize=7.5,
           color='0.3', ha='right')
ax[1].set_xlabel('RMS timing jitter [chip periods]')
ax[1].set_ylabel('RMS Bragg error [pm]')
ax[1].set_title('(b) Timing jitter tolerance')
ax[1].legend(fontsize=7.5); ax[1].grid(True, which='both', alpha=0.25)
plt.tight_layout(); plt.savefig('figs/fig_s13_chiprate.png', dpi=140)

# ---------------------------------------------------------------------------
# printed tables
# ---------------------------------------------------------------------------
print('acquisition-chain study: K=%d gratings, N=%d chips, M=%d slots, R=%.2f'
      % (K, N_CHIPS, M_SLOTS, R))
print('--- (a) RMS Bragg error [pm] vs chip rate ---')
print('   B[Mchip/s]   dz[m]   d_min=0.5 m   d_min=2 m   d_min=4 m')
for i, B in enumerate(Bs):
    print('   %10.0f   %5.2f   %11.2f   %9.2f   %9.2f'
          % (B / 1e6, dz_of(B), e_05[i], e_2[i], e_4[i]))
for d, e in [(0.5, e_05), (2.0, e_2), (4.0, e_4)]:
    ok = Bs[e <= 10.0]
    req = C_LIGHT / (2 * N_GROUP * d) / 1e6
    print('   d_min = %.1f m: criterion dz = d_min gives B = %.1f Mchip/s;'
          ' simulation meets 10 pm from B = %s Mchip/s'
          % (d, req, ('%.0f' % (ok[0] / 1e6)) if len(ok) else 'never'))
print('--- (b) timing jitter tolerance at B = %.0f Mchip/s ---' % (B_WORK / 1e6))
print('   (array with 4 m minimum spacing)')
print('   jitter[chip]   jitter@25M[ns]   err B=25M [pm]   err B=50M [pm]')
for j, a, b in zip(jit, e_jit_25, e_jit_50):
    print('   %12.2f   %14.2f   %14.2f   %14.2f' % (j, j / 25e6 * 1e9, a, b))
print('--- equivalent-time sampling budget (OSR = %d samples/chip) ---' % OSR)
print('   B[Mchip/s]   ADC[MS/s]   samples/period   periods   T_acq[us]   step[ns]')
for B in [25e6, 50e6, 100e6]:
    for fs in [3.6e6, 20e6]:
        t_acq, per, periods, step = ets(B, fs)
        print('   %10.0f   %9.1f   %14d   %7d   %9.1f   %8.2f'
              % (B / 1e6, fs / 1e6, per, periods, t_acq * 1e6, step * 1e9))
print('--- converter resolution budget ---')
print('   the grating peak rarely fills the converter: the unipolar drive and the')
print('   returns of all the other gratings sit under it as a pedestal')
print('   bits   jitter @50%% headroom [pm]   jitter @2%% headroom [pm]')
for b in [8, 10, 12, 14, 16]:
    j50, _ = quant_jitter(b, headroom=0.5)
    j02, _ = quant_jitter(b, headroom=0.02)
    print('   %4d   %24.3f   %23.3f' % (b, j50, j02))
print('saved figs/fig_s13_chiprate.png')
