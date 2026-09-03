"""s52_schemes.py - code against chirp, pulse and CW interrogation, numerically.

Four ways of reading the same serial array with the same swept laser, the same
peak optical power and the same measurement time per wavelength step (one
period of N chips):

  CDM    unipolar m-sequence of N chips, correlation with the bipolar replica
         after mean removal (the receiver of the paper), constant offset of
         the periodic correlation subtracted.
  chirp  intensity chirp 0.5 (1 + cos phi(t)) with the instantaneous frequency
         swept linearly over 0.02 B .. B in one period, matched filter with a
         Hamming-weighted replica (rectangular weighting would leave -13 dB
         side lobes).
  pulse  one chip of light per period, matched filter with the chip (a boxcar).
  CW     constant power, no delay information: the record mean is the summed
         reflection of every grating in the band, so only K = 1 per band is
         readable. Its error is reported for K = 1 (the ideal case) and, for
         reference, for K = 3 read as if the summed spectrum were each
         grating's own.

Every scheme runs through the time-domain model of s50: three gratings at 8,
19.2 and 40 m, every path up to the third order with the (1-R) transmission of
every grating crossed, spectral shadowing through the wavelength dependence of
R, link-budget noise (P0 into 4 dB of loss, NEP 0.5 pW/sqrt(Hz), shot noise,
RIN -130 dB/Hz), a 4th-order Bessel low-pass at 0.75 B, a 12-bit ADC at 8
samples per chip. The wavelength sweep, the detunings and the peak estimator
are those of Table III (M = 64 steps over +-2.6 FWHM, detunings within +-200
pm, Gaussian fit to the samples above 30 % within +-14 % of the sweep).

Output: RMS wavelength error per scheme and case, and the SNR of the delay
peak of the first grating at its line centre.
"""
import os
import numpy as np
from scipy.signal import bessel, filtfilt
import common as C

C_LIGHT, NG = 2.998e8, 1.468
CHIP_RATE, NBITS = 25e6, 7
SPC_FINE, NS = 32, 8
Z = {'b': 8.0, 'c': 19.2, 'a': 40.0}
ORDER = sorted(Z, key=Z.get)
M_PER_CHIP = C_LIGHT / NG / CHIP_RATE / 2.0
ALPHA, NEP, RESP, RIN_DB, ADC_BITS, Q_E = 10 ** (-4.0 / 10), 0.5e-12, 0.9, -130.0, 12, 1.602e-19
FWHM_GHZ = C.FBG_FWHM_GHZ
PM = 1.0 / C.GHZ_PER_PM                       # pm per GHz
M_STEPS = 64
NU = np.linspace(-2.6 * FWHM_GHZ, 2.6 * FWHM_GHZ, M_STEPS)
DETUNE_PM = 200.0
N_DRAWS = 5

CODE01 = C._mls01(NBITS)
N = CODE01.size
FS = CHIP_RATE * SPC_FINE
T_FINE = np.arange(N * SPC_FINE) / FS
B_LPF, A_LPF = bessel(4, 0.75 * CHIP_RATE / (FS / 2), norm='mag')


def line(nu, nu0):
    return np.exp(-0.5 * ((nu - nu0) / (FWHM_GHZ / 2.35482)) ** 2)


# ---------------------------------------------------------------------------
# transmitted waveforms (fraction of the peak power) and matched replicas
# ---------------------------------------------------------------------------
def waveform(scheme):
    if scheme == 'CDM':
        tx = np.repeat(CODE01.astype(float), SPC_FINE)
        rep_ = np.repeat(C._to_pm1(CODE01).astype(float), SPC_FINE)
        return tx, rep_
    if scheme == 'chirp':
        T = N / CHIP_RATE
        f0, f1 = 0.02 * CHIP_RATE, CHIP_RATE
        phi = 2 * np.pi * (f0 * T_FINE + 0.5 * (f1 - f0) / T * T_FINE ** 2)
        tx = 0.5 * (1.0 + np.cos(phi))
        rep_ = np.cos(phi) * np.hamming(T_FINE.size)
        return tx, rep_
    if scheme == 'pulse':
        tx = np.zeros(N * SPC_FINE); tx[:SPC_FINE] = 1.0
        return tx, tx.copy()
    if scheme == 'CW':
        tx = np.ones(N * SPC_FINE)
        return tx, None
    raise ValueError(scheme)


# ---------------------------------------------------------------------------
# every path up to the third order at one wavelength
# ---------------------------------------------------------------------------
def paths(Rl):
    """Rl: reflectivity of every grating at the current wavelength."""
    def trans(z0, z1):
        t = 1.0
        for g in ORDER:
            if min(z0, z1) < Z[g] < max(z0, z1):
                t *= (1.0 - Rl[g])
        return t

    def one(seq):
        amp, zz, length = 1.0, 0.0, 0.0
        for g in seq:
            amp *= trans(zz, Z[g]) * Rl[g]
            length += abs(Z[g] - zz)
            zz = Z[g]
        amp *= trans(zz, 0.0)
        return amp, (length + zz) / 2.0

    out = [((g,),) + one((g,)) for g in ORDER]
    for g1 in ORDER:
        for g2 in ORDER:
            for g3 in ORDER:
                if Z[g2] < Z[g1] and Z[g2] < Z[g3]:
                    out.append(((g1, g2, g3),) + one((g1, g2, g3)))
    return out


def record(scheme, Rl, P0, rng, noise=True, gratings=ORDER):
    """ADC samples of one period at one wavelength, in amperes."""
    tx, _ = waveform(scheme)
    popt = np.zeros_like(tx)
    for seq, amp, zpos in paths(Rl):
        if any(g not in gratings for g in seq):
            continue
        popt += amp * np.roll(tx, int(round(zpos / M_PER_CHIP * SPC_FINE)))
    popt *= P0 * ALPHA
    if noise:
        popt = popt * (1.0 + rng.normal(0.0, np.sqrt(10 ** (RIN_DB / 10) * FS / 2), popt.size))
    i_pd = RESP * popt
    if noise:
        i_pd = i_pd + rng.normal(0.0, np.sqrt(2 * Q_E * RESP * max(popt.mean(), 0.0) * FS / 2), i_pd.size)
        i_pd = i_pd + rng.normal(0.0, RESP * NEP * np.sqrt(FS / 2), i_pd.size)
    i_pd = filtfilt(B_LPF, A_LPF, i_pd)
    dec = SPC_FINE // NS
    return i_pd[dec // 2::dec]


def quantize(rec, fsr):
    return np.round(rec / fsr * 2 ** (ADC_BITS - 1)) / 2 ** (ADC_BITS - 1) * fsr


def estimate(scheme, rec):
    """Amplitude of the delay profile at every grating (a.u.), or the record
    mean for CW."""
    if scheme == 'CW':
        return {g: rec.mean() for g in ORDER}
    _, rep_ = waveform(scheme)
    dec = SPC_FINE // NS
    rep_ = rep_[dec // 2::dec]
    r = rec - rec.mean()
    corr = np.fft.ifft(np.fft.fft(r) * np.conj(np.fft.fft(rep_))).real
    if scheme == 'CDM':
        corr = -corr                              # unipolar code against the bipolar replica
    corr -= np.median(corr)                       # constant offset of the profile (the mean removal leaks -1/N of every return into every bin)
    zaxis = np.arange(rec.size) / NS * M_PER_CHIP
    out = {}
    half = 0.6 * M_PER_CHIP
    for g in ORDER:
        w = np.abs(zaxis - Z[g]) <= half
        out[g] = corr[w].max()
    return out


def sweep(scheme, R, P0, det_pm, rng, gratings=ORDER, noise=True):
    """Spectrum of every grating (M_STEPS values) as the scheme sees it."""
    nu0 = {g: det_pm[g] / PM for g in ORDER}
    # full scale of the ADC from a noise-free record at the strongest step
    Rl_max = {g: R for g in ORDER}
    fsr = 1.2 * np.abs(record(scheme, Rl_max, P0, rng, noise=False, gratings=gratings)).max()
    spec = {g: np.zeros(M_STEPS) for g in ORDER}
    for i, nu in enumerate(NU):
        Rl = {g: R * line(nu, nu0[g]) for g in ORDER}
        rec = quantize(record(scheme, Rl, P0, rng, noise=noise, gratings=gratings), fsr)
        est = estimate(scheme, rec)
        for g in ORDER:
            spec[g][i] = est[g]
    return spec, nu0


def rms_error(scheme, R, P0, rng, gratings=ORDER):
    errs = []
    for d in range(N_DRAWS):
        det = {g: rng.uniform(-DETUNE_PM, DETUNE_PM) for g in ORDER}
        spec, nu0 = sweep(scheme, R, P0, det, rng, gratings=gratings)
        for g in gratings:
            mu = C.gauss_fit_peak(NU, spec[g]) * PM
            errs.append(mu - det[g])
    return float(np.sqrt(np.mean(np.square(errs))))


def snr_db(scheme, R, P0, rng):
    """Peak of the first grating at its line centre against the noise there."""
    Rl = {g: R for g in ORDER}
    rec0 = record(scheme, Rl, P0, rng, noise=False)
    fsr = 1.2 * np.abs(rec0).max()
    e0 = estimate(scheme, quantize(rec0, fsr))
    vals = []
    for _ in range(20):
        e = estimate(scheme, quantize(record(scheme, Rl, P0, rng, noise=True), fsr))
        vals.append(e[ORDER[0]])
    return 20 * np.log10(e0[ORDER[0]] / max(np.std(vals), 1e-30))


if __name__ == '__main__':
    CASES = [(0.0, 0.10), (0.0, 0.01), (0.0, 0.001), (-20.0, 0.001)]
    SCHEMES = ['CDM', 'chirp', 'pulse', 'CW']
    rng = np.random.default_rng(11)
    lines = []
    hdr = 'P0 [dBm]  R      ' + ''.join('%14s' % s for s in SCHEMES) + '   CW, K=1'
    lines.append('RMS wavelength error [pm], K = 3 gratings (CW: summed spectrum read as each grating)')
    lines.append(hdr)
    for p0_dbm, R in CASES:
        P0 = 10 ** (p0_dbm / 10) * 1e-3
        row = '%8.0f  %5.1f%% ' % (p0_dbm, 100 * R)
        for s in SCHEMES:
            row += '%14.2f' % rms_error(s, R, P0, np.random.default_rng(11))
        row += '%10.2f' % rms_error('CW', R, P0, np.random.default_rng(11), gratings=(ORDER[0],))
        lines.append(row)
        print(row, flush=True)
    lines.append('')
    lines.append('SNR of the delay peak of grating b at its line centre [dB, 20 log]')
    lines.append(hdr)
    for p0_dbm, R in CASES:
        P0 = 10 ** (p0_dbm / 10) * 1e-3
        row = '%8.0f  %5.1f%% ' % (p0_dbm, 100 * R)
        for s in SCHEMES:
            row += '%14.1f' % snr_db(s, R, P0, np.random.default_rng(5))
        lines.append(row)
        print(row, flush=True)
    os.makedirs('out', exist_ok=True)
    with open('out/s52_schemes.txt', 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
