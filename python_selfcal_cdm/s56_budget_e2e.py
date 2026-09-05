"""s56_budget_e2e.py - one consistent run of the whole readout, and the error
budget derived from it by switching mechanisms off one at a time.

Why. The earlier budget (s15) added rows that came from separate models, and
its 'delay-bin resolution' row turned out to be code leakage read in a rounded
bin. Here every number comes from the same chain, for the same random draws:

  optics    K gratings on one fiber, Gaussian lines with a random edge
            asymmetry, random detunings within +-W, spectral shadowing by the
            product of transmissions, every third-order ghost path with its own
            transmissions and its delay wrapped modulo the code period,
            references on a separate stub (no shadowing between them and the
            sensors).
  source    stepped sweep of M steps. At each step the laser dwells on a
            distribution p(xi) of wavelength offsets (skew-normal, standard
            deviation Delta_ch/2, mean offset that varies smoothly across the
            band). The whole received spectrum is averaged over p, sensors and
            ghosts alike, so chirp and shadowing are not separated by hand.
  receiver  one code period of the unipolar m-sequence at SPC_FINE samples per
            chip, link loss, laser RIN, shot noise, NEP, 4th-order Bessel
            low-pass at 0.75 B, NS ADC samples per chip, 12-bit quantization,
            mean removal, periodic correlation with the bipolar replica. Two
            receivers: 'baseline' (leaves the -1/N side lobes) and 'offset'
            (subtracts the median of the profile, Sec. III-D).
  readout   the correlation is read at the true delay (linear interpolation)
            or at the nearest chip bin ('bin' readout). Gaussian LS fit with
            the estimator of Table III. Optional sequential deshadowing with
            fitted transmissions and a 0.05 floor, as in s12/s15. The
            wavelength axis carries a table drift (offset + gain) that the
            references remove with a polynomial of order nref-1. Reference
            readings get a 1-pm RMS stability error.
  truth     the fitted center of the same grating read alone, without chirp,
            through the ideal chain. That removes the static line-shape
            offset, which the sensor calibration removes in practice.

Ablations ('without X') switch one mechanism off and keep everything else,
including the random draws. The difference to the full run is the price of
that mechanism in this configuration. 'Only X' keeps one mechanism and
switches the others off.

Output: out/s56_budget_e2e.txt and out/s56_budget_e2e.npz.
Runtime: a few minutes per configuration at NTRIALS = 20.
"""
import os, sys, time, warnings
import numpy as np
from scipy.signal import bessel, filtfilt
from scipy.special import erf
import common as C
warnings.filterwarnings('ignore')

PM = C.PM_PER_GHZ
C_LIGHT, NG = 2.998e8, 1.468
SPC_FINE, NS, ADC_BITS = 16, 8, 12
P0, ALPHA, NEP, RESP, RIN_DB = 1e-3, 10 ** (-4.0 / 10), 0.5e-12, 0.9, -130.0
Q_E = 1.602e-19
W_GHZ = 200.0 / PM                    # half band, +-200 pm
M_STEPS = 64
LAM = np.linspace(-650.0, 650.0, M_STEPS) / PM     # sweep steps, GHz
XI = np.linspace(-4.0, 4.0, 33)       # chirp kernel grid in units of its std
REF_STAB_PM = 1.0
DRIFT_OFF_PM, DRIFT_GAIN = 20.0, 0.003
ASYM_SPREAD = 0.30
NTRIALS = int(os.environ.get('S56_NTRIALS', '20'))

CFG = {
    'bench':    dict(K=3, R=0.10, fwhm=250.0, nbits=7, B=25e6, spacing='uniform', d_min=4.0,
                     nref=2, peel=False, chirp=0.2),
    'designed': dict(K=32, R=0.01, fwhm=250.0, nbits=9, B=100e6, spacing='random', d_min=2.0,
                     nref=3, peel=True, chirp=0.2),
}

ALL = ('shadow', 'ghosts', 'leak', 'noise', 'chirp', 'resolution', 'refs', 'drift')


def kernel(std_ghz, skew=1.2):
    if std_ghz <= 0:
        return np.array([0.0]), np.array([1.0])
    d = XI * std_ghz
    p = np.exp(-0.5 * XI ** 2) * (1.0 + erf(skew * XI / np.sqrt(2)))
    p = np.clip(p, 0, None)
    p /= p.sum()
    return d, p


def make_layout(cfg, rng):
    K, N = cfg['K'], 2 ** cfg['nbits'] - 1
    m_per_chip = C_LIGHT / NG / cfg['B'] / 2.0
    if cfg['spacing'] == 'uniform':
        z = cfg['d_min'] * np.arange(1, K + 1)
    elif cfg['spacing'] == 'given':
        z = np.asarray(cfg['z'], float)
    else:
        # random bins, at least d_min apart, inside 80 percent of the period
        while True:
            zc = np.sort(rng.uniform(2.0, 0.8 * N * m_per_chip, K))
            if np.all(np.diff(zc) >= cfg['d_min']):
                break
        z = zc
    tau = z / m_per_chip                      # chips, fractional
    # references on a stub: far bins, at least 3 chips from every sensor and each other
    refs = []
    cand = np.arange(2, N - 2)
    while len(refs) < cfg['nref']:
        b = rng.choice(cand)
        if np.all(np.abs(b - tau) >= 3.0) and all(abs(b - r) >= 3.0 for r in refs):
            refs.append(float(b))
    return tau, np.array(refs)


def spectra_matrix(cfg, nub, asym, tau_s, tau_r, nub_r, on, std_chirp, off_fn):
    """Return (delay list [chips], amplitude matrix [paths x M]) for all paths.
    Sensors: direct returns with shadowing and every third-order ghost. The
    references sit on a stub, so they only see themselves."""
    K = len(nub)
    F = cfg['fwhm'] / PM
    R = cfg['R']
    dxi, pxi = kernel(std_chirp) if on['chirp'] else (np.array([0.0]), np.array([1.0]))
    # evaluation grid: steps x kernel offsets (chirp mean offset varies with the band position of ... the step)
    lam = LAM[:, None] + dxi[None, :]                      # M x X
    if on['chirp']:
        lam = lam + off_fn(LAM)[:, None]
    shape = lambda nb, a: C.fbg_tanh(lam, nb, F, n_side=a) if False else \
        np.exp(-0.5 * ((lam - nb) / (F / 2.35482)) ** 2) * (1.0 + a * np.tanh((lam - nb) / (F / 2.35482)))
    Rk = np.array([np.clip(R * shape(nub[k], asym[k]), 0, None) for k in range(K)])   # K x M x X
    Rk = Rk / np.max(Rk, axis=(1, 2), keepdims=True) * R
    if on['shadow']:
        logT = np.log(np.clip(1.0 - Rk, 1e-6, None))       # one-way log transmission
    else:
        logT = np.zeros_like(Rk)
    cumL = np.concatenate([np.zeros((1,) + Rk.shape[1:]), np.cumsum(logT, axis=0)], axis=0)  # L_x = sum over index < x
    delays, amps = [], []
    for k in range(K):
        amps.append(Rk[k] * np.exp(2.0 * cumL[k]))          # out and back through gratings < k
        delays.append(tau_s[k])
    if on['ghosts']:
        idx = np.arange(K)
        a, b, c = np.meshgrid(idx, idx, idx, indexing='ij')
        m = (b < a) & (b < c)
        a, b, c = a[m], b[m], c[m]
        for ai, bi, ci in zip(a, b, c):
            # 0 -> a -> b -> c -> 0, transmissions strictly between the endpoints
            T = np.exp(cumL[ai] + (cumL[ai] - cumL[bi + 1]) + (cumL[ci] - cumL[bi + 1]) + cumL[ci])
            amps.append(Rk[ai] * Rk[bi] * Rk[ci] * T)
            delays.append(tau_s[ai] - tau_s[bi] + tau_s[ci])
    # references: alone on their stub
    for j in range(len(tau_r)):
        amps.append(np.clip(R * shape(nub_r[j], 0.0), 0, None))
        delays.append(tau_r[j])
    A = np.array([(am * pxi[None, :]).sum(axis=1) for am in amps])    # paths x M, chirp-averaged
    return np.array(delays), A


def receiver(cfg, delays, A, on, rng, mode):
    """Correlation profiles for every step: returns profile [M x N*NS] in units
    where a unit return gives a peak of 1."""
    N = 2 ** cfg['nbits'] - 1
    code01 = C._mls01(cfg['nbits']).astype(float)
    B = cfg['B']
    fs = B * SPC_FINE
    L = N * SPC_FINE
    tx = np.repeat(code01, SPC_FINE)
    TX = np.fft.fft(tx)
    replica = np.repeat(C._to_pm1(code01).astype(float), NS)
    REP = np.conj(np.fft.fft(replica))
    b_, a_ = bessel(4, 0.75 * B / (fs / 2), norm='mag')
    pos = np.mod(np.round(delays * SPC_FINE).astype(int), L)
    out = np.zeros((A.shape[1], N * NS))
    for m in range(A.shape[1]):
        hist = np.zeros(L)
        np.add.at(hist, pos, A[:, m])
        popt = np.fft.ifft(np.fft.fft(hist) * TX).real * P0 * ALPHA
        if on['noise']:
            popt = popt * (1.0 + rng.normal(0.0, np.sqrt(10 ** (RIN_DB / 10) * fs / 2), L))
        i_pd = RESP * popt
        if on['noise']:
            i_pd = i_pd + rng.normal(0.0, np.sqrt(2 * Q_E * RESP * max(popt.mean(), 0) * fs / 2), L)
            i_pd = i_pd + rng.normal(0.0, RESP * NEP * np.sqrt(fs / 2), L)
        if on['leak']:
            i_pd = filtfilt(b_, a_, i_pd)
            dec = SPC_FINE // NS
            rec = i_pd[dec // 2::dec]
            fsr = 1.2 * max(rec.max(), 1e-12)
            rec = np.round(rec / fsr * 2 ** (ADC_BITS - 1)) / 2 ** (ADC_BITS - 1) * fsr
            rec = rec - rec.mean()
            corr = np.fft.ifft(np.fft.fft(rec) * REP).real
            corr /= -(N * NS / 2.0) * RESP * P0 * ALPHA
            if mode == 'offset':
                corr -= np.median(corr)
        else:
            # ideal correlation: the delay histogram itself, seen through the same
            # chip shape (triangular overlap), without side lobes, noise-free chain
            dec = SPC_FINE // NS
            h = np.zeros(N * NS)
            np.add.at(h, np.mod(np.round(delays * NS).astype(int), N * NS), A[:, m])
            tri = np.concatenate([np.arange(NS, 0, -1), np.zeros(N * NS - 2 * NS + 1), np.arange(1, NS)]) / NS
            corr = np.fft.ifft(np.fft.fft(h) * np.fft.fft(tri)).real
        out[m] = corr
    return out


def read(profile, tau_chips, resolution_on):
    """Sample S(lambda_m) at a delay: nearest chip bin ('resolution' on) or the true delay."""
    n = profile.shape[1]
    if resolution_on:
        x = np.round(tau_chips) * NS
    else:
        x = tau_chips * NS
    i0 = int(np.floor(x)) % n
    f = x - np.floor(x)
    return (1 - f) * profile[:, i0] + f * profile[:, (i0 + 1) % n]


def peel(S, nub_guess, floor=0.05):
    """sequential deshadowing with fitted Gaussian transmissions, as in s12/s15"""
    K = S.shape[0]
    T = np.ones(M_STEPS)
    cent = np.empty(K)
    for k in range(K):
        y = S[k] / np.maximum(T, floor)
        amp, mu, sg, _ = C.gauss_fit_full(LAM, y)
        cent[k] = mu
        T = T * (1.0 - np.clip(amp, 0, 0.99) * np.exp(-0.5 * ((LAM - mu) / max(sg, 1e-3)) ** 2)) ** 2
    return cent


def one_trial(cfg, seed, on, mode):
    rng = np.random.default_rng(seed)
    K = cfg['K']
    tau_s, tau_r = make_layout(cfg, rng)
    if 'nub_pm' in cfg:
        nub = np.asarray(cfg['nub_pm'], float) / PM
    else:
        nub = rng.uniform(-W_GHZ, W_GHZ, K)
    asym = rng.uniform(-ASYM_SPREAD, ASYM_SPREAD, K)
    nub_r = np.linspace(-0.9 * W_GHZ, 0.9 * W_GHZ, cfg['nref']) if cfg['nref'] > 1 else np.array([0.0])
    F = cfg['fwhm'] / PM
    std_chirp = cfg['chirp'] * F / 2.0
    off_fn = lambda nb: std_chirp * (0.20 + 0.20 * nb / W_GHZ + 0.15 * (nb / W_GHZ) ** 2 + 0.10 * np.sin(2.5 * nb / W_GHZ))
    delays, A = spectra_matrix(cfg, nub, asym, tau_s, tau_r, nub_r, on, std_chirp, off_fn)
    prof = receiver(cfg, delays, A, on, rng, mode)
    S = np.array([read(prof, t, on['resolution']) for t in tau_s])
    Sr = np.array([read(prof, t, on['resolution']) for t in tau_r])
    # axis: table drift, applied to every reading (the laser is where the table says plus drift)
    drift = (DRIFT_OFF_PM + DRIFT_GAIN * LAM * PM) / PM if on['drift'] else np.zeros(M_STEPS)
    lam_axis = LAM - drift        # the sweep visits LAM, but the table reports LAM - drift ... handled by fitting on reported axis
    if cfg['peel'] and on['shadow']:
        cent = peel(S, nub)
    else:
        cent = np.array([C.gauss_fit_peak(LAM, S[k]) for k in range(K)])
    cent_r = np.array([C.gauss_fit_peak(LAM, Sr[j]) for j in range(cfg['nref'])])
    # reported values on the table axis
    cent_rep = cent - np.interp(cent, LAM, drift)
    cent_r_rep = cent_r - np.interp(cent_r, LAM, drift)
    # reference correction: known nub_r, readings carry a stability error
    if on['refs']:
        stab = rng.normal(0.0, REF_STAB_PM / PM, cfg['nref'])
        err_r = cent_r_rep - (nub_r + stab)
        p = np.polyfit(nub_r, err_r, cfg['nref'] - 1) if cfg['nref'] > 1 else np.array([err_r.mean()])
        cent_rep = cent_rep - np.polyval(p, cent_rep)
    # truth: each sensor alone, no chirp, ideal chain, static line-shape offset removed
    truth = np.empty(K)
    for k in range(K):
        y = np.exp(-0.5 * ((LAM - nub[k]) / (F / 2.35482)) ** 2) * (1.0 + asym[k] * np.tanh((LAM - nub[k]) / (F / 2.35482)))
        truth[k] = C.gauss_fit_peak(LAM, np.clip(y, 0, None))
    return (cent_rep - truth) * PM


def run(cfg, on, mode, ntrials):
    errs = np.concatenate([one_trial(cfg, 1000 + t, on, mode) for t in range(ntrials)])
    return errs


def stats(e):
    a = np.abs(e)
    return dict(rms=float(np.sqrt(np.mean(e ** 2))), p50=float(np.median(a)), p90=float(np.percentile(a, 90)), mx=float(a.max()))


if __name__ == '__main__':
    lines = []
    say = lambda s='': (print(s), lines.append(s))
    results = {}
    t0 = time.time()
    for name, cfg in CFG.items():
        for mode in ('baseline', 'offset'):
            full_on = {k: True for k in ALL}
            e_full = run(cfg, full_on, mode, NTRIALS)
            st = stats(e_full)
            say('=== %s, %s receiver, %d trials x %d sensors: RMS %.2f pm, |e| p50 %.2f, p90 %.2f, max %.2f'
                % (name, mode, NTRIALS, cfg['K'], st['rms'], st['p50'], st['p90'], st['mx']))
            results['%s/%s/full' % (name, mode)] = e_full
            say('%-14s %10s %10s' % ('mechanism', 'without', 'only'))
            for mech in ALL:
                on_w = dict(full_on); on_w[mech] = False
                e_w = run(cfg, on_w, mode, NTRIALS)
                on_o = {k: (k == mech) for k in ALL}
                # 'only' keeps the reference correction and the drift together (they are a pair)
                if mech in ('refs', 'drift'):
                    on_o['refs'] = True; on_o['drift'] = True
                e_o = run(cfg, on_o, mode, NTRIALS)
                results['%s/%s/without_%s' % (name, mode, mech)] = e_w
                results['%s/%s/only_%s' % (name, mode, mech)] = e_o
                say('%-14s %10.2f %10.2f' % (mech, stats(e_w)['rms'], stats(e_o)['rms']))
            say('   (%.0f s)' % (time.time() - t0))
            say()
    os.makedirs('out', exist_ok=True)
    open('out/s56_budget_e2e.txt', 'w').write('\n'.join(lines) + '\n')
    np.savez('out/s56_budget_e2e.npz', **{k.replace('/', '__'): v for k, v in results.items()})
    print('saved out/s56_budget_e2e.txt and .npz')
