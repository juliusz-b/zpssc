"""s58_decorrelator.py - a decorrelating receiver instead of the matched filter.

The correlation receiver reads grating k at its delay and accepts whatever the
code side lobes and the neighbouring correlation peaks put there. But the
mixing is linear and known: the record is y = D a + n, where the columns of D
are the delayed code waveforms of the known gratings after the receiver
filter and ADC, and a are the amplitudes A_k(lambda_m). Solving the least
squares problem a = (D'D)^-1 D' y removes the side lobes and the overlap of
neighbouring peaks exactly, for any code, at a noise penalty set by the
conditioning of D. This is the zero-forcing multiuser detector of CDMA,
applied to the FBG array.

Tested on the full chain of s56 (same optics, chirp, noise, references): the
bench (three 10 percent gratings 0.98 chip apart, where the overlap of
correlation peaks leaves 2.5 pm even in an ideal chain) and the designed
array. Then the leakage-limited regime: K = 32 .. 120 gratings at N = 127,
R = 0.3 percent, where the correlation receiver is capped at K ~ 47 by the
-1/N side lobes (Fig. 11), against the decorrelator.

Output: out/s58_decorrelator.txt
"""
import os, time, warnings
import numpy as np
from scipy.signal import bessel, filtfilt
import common as C
import s56_budget_e2e as S
warnings.filterwarnings('ignore')

NT = int(os.environ.get('S58_NTRIALS', '20'))


def design_matrix(cfg, delays_known):
    """columns: unit returns at the known delays, through the same filter and ADC"""
    N = 2 ** cfg['nbits'] - 1
    code01 = C._mls01(cfg['nbits']).astype(float)
    B = cfg['B']; fs = B * S.SPC_FINE; L = N * S.SPC_FINE
    tx = np.repeat(code01, S.SPC_FINE); TX = np.fft.fft(tx)
    b_, a_ = bessel(4, 0.75 * B / (fs / 2), norm='mag')
    dec = S.SPC_FINE // S.NS
    cols = []
    for d in delays_known:
        h = np.zeros(L); h[int(np.round(d * S.SPC_FINE)) % L] = 1.0
        w = np.fft.ifft(np.fft.fft(h) * TX).real
        w = filtfilt(b_, a_, w)[dec // 2::dec]
        cols.append(w - w.mean())          # the receiver removes the record mean
    return np.array(cols).T                # (N*NS) x K


def receiver_lsq(cfg, delays, A, on, rng, delays_known):
    """same front end as s56.receiver, then least squares instead of correlation"""
    N = 2 ** cfg['nbits'] - 1
    code01 = C._mls01(cfg['nbits']).astype(float)
    B = cfg['B']; fs = B * S.SPC_FINE; L = N * S.SPC_FINE
    tx = np.repeat(code01, S.SPC_FINE); TX = np.fft.fft(tx)
    b_, a_ = bessel(4, 0.75 * B / (fs / 2), norm='mag')
    dec = S.SPC_FINE // S.NS
    D = design_matrix(cfg, delays_known)
    pinv = np.linalg.pinv(D)
    pos = np.mod(np.round(delays * S.SPC_FINE).astype(int), L)
    est = np.zeros((A.shape[1], len(delays_known)))
    for m in range(A.shape[1]):
        hist = np.zeros(L); np.add.at(hist, pos, A[:, m])
        popt = np.fft.ifft(np.fft.fft(hist) * TX).real * S.P0 * S.ALPHA
        if on['noise']:
            popt = popt * (1.0 + rng.normal(0.0, np.sqrt(10 ** (S.RIN_DB / 10) * fs / 2), L))
        i_pd = S.RESP * popt
        if on['noise']:
            i_pd = i_pd + rng.normal(0.0, np.sqrt(2 * S.Q_E * S.RESP * max(popt.mean(), 0) * fs / 2), L)
            i_pd = i_pd + rng.normal(0.0, S.RESP * S.NEP * np.sqrt(fs / 2), L)
        i_pd = filtfilt(b_, a_, i_pd)
        rec = i_pd[dec // 2::dec]
        fsr = 1.2 * max(rec.max(), 1e-12)
        rec = np.round(rec / fsr * 2 ** (S.ADC_BITS - 1)) / 2 ** (S.ADC_BITS - 1) * fsr
        rec = rec - rec.mean()
        est[m] = pinv @ rec / (S.RESP * S.P0 * S.ALPHA)
    return est                                   # M x K_known


def one_trial_lsq(cfg, seed, on):
    """s56.one_trial with the decorrelating receiver, references included as known delays"""
    rng = np.random.default_rng(seed)
    K = cfg['K']
    tau_s, tau_r = S.make_layout(cfg, rng)
    nub = np.asarray(cfg['nub_pm'], float) / S.PM if 'nub_pm' in cfg else rng.uniform(-S.W_GHZ, S.W_GHZ, K)
    asym = rng.uniform(-S.ASYM_SPREAD, S.ASYM_SPREAD, K)
    nub_r = np.linspace(-0.9 * S.W_GHZ, 0.9 * S.W_GHZ, cfg['nref']) if cfg['nref'] > 1 else np.array([0.0])
    F = cfg['fwhm'] / S.PM
    std_chirp = cfg['chirp'] * F / 2.0
    off_fn = lambda nb: std_chirp * (0.20 + 0.20 * nb / S.W_GHZ + 0.15 * (nb / S.W_GHZ) ** 2 + 0.10 * np.sin(2.5 * nb / S.W_GHZ))
    delays, A = S.spectra_matrix(cfg, nub, asym, tau_s, tau_r, nub_r, on, std_chirp, off_fn)
    known = np.concatenate([tau_s, tau_r])
    est = receiver_lsq(cfg, delays, A, on, rng, known)
    Ssen = est[:, :K].T; Sref = est[:, K:].T
    drift = (S.DRIFT_OFF_PM + S.DRIFT_GAIN * S.LAM * S.PM) / S.PM if on['drift'] else np.zeros(S.M_STEPS)
    if cfg['peel'] and on['shadow']:
        cent = S.peel(Ssen, nub)
    else:
        cent = np.array([C.gauss_fit_peak(S.LAM, Ssen[k]) for k in range(K)])
    cent_r = np.array([C.gauss_fit_peak(S.LAM, Sref[j]) for j in range(cfg['nref'])])
    cent_rep = cent - np.interp(cent, S.LAM, drift)
    cent_r_rep = cent_r - np.interp(cent_r, S.LAM, drift)
    if on['refs']:
        stab = rng.normal(0.0, S.REF_STAB_PM / S.PM, cfg['nref'])
        err_r = cent_r_rep - (nub_r + stab)
        p = np.polyfit(nub_r, err_r, cfg['nref'] - 1) if cfg['nref'] > 1 else np.array([err_r.mean()])
        cent_rep = cent_rep - np.polyval(p, cent_rep)
    truth = np.empty(K)
    for k in range(K):
        y = np.exp(-0.5 * ((S.LAM - nub[k]) / (F / 2.35482)) ** 2) * (1.0 + asym[k] * np.tanh((S.LAM - nub[k]) / (F / 2.35482)))
        truth[k] = C.gauss_fit_peak(S.LAM, np.clip(y, 0, None))
    return (cent_rep - truth) * S.PM


if __name__ == '__main__':
    lines = []
    say = lambda s='': (print(s), lines.append(s))
    on = {k: True for k in S.ALL}
    t0 = time.time()
    say('=== full chain, 20 trials: matched filter (baseline), offset subtraction, decorrelator')
    for name in ('bench', 'designed'):
        cfg = S.CFG[name]
        e_b = S.run(cfg, on, 'baseline', NT); e_o = S.run(cfg, on, 'offset', NT)
        e_l = np.concatenate([one_trial_lsq(cfg, 1000 + t, on) for t in range(NT)])
        sb, so, sl = S.stats(e_b), S.stats(e_o), S.stats(e_l)
        say('%-9s RMS  baseline %5.2f  offset %5.2f  decorrelator %5.2f pm   | p90 %5.2f %5.2f %5.2f | max %5.2f %5.2f %5.2f'
            % (name, sb['rms'], so['rms'], sl['rms'], sb['p90'], so['p90'], sl['p90'], sb['mx'], so['mx'], sl['mx']))
        on_ns = dict(on); on_ns['noise'] = False
        e_ln = np.concatenate([one_trial_lsq(cfg, 1000 + t, on_ns) for t in range(NT)])
        say('%-9s decorrelator without noise %.2f pm (the noise penalty is the difference)' % (name, S.stats(e_ln)['rms']))
    say('   (%.0f s)' % (time.time() - t0))
    say()
    say('=== leakage-limited regime: N = 127, R = 0.3 percent, distinct chip bins plus a random fraction, no chirp, 8 trials')
    say('%6s %12s %12s %12s' % ('K', 'baseline', 'offset', 'decorrelator'))
    for K in (16, 32, 48, 64, 80, 96, 110):
        cfg = dict(K=K, R=0.003, fwhm=250.0, nbits=7, B=100e6, spacing='bins', d_min=1.0,
                   nref=2, peel=False, chirp=0.0)
        try:
            e_b = S.run(cfg, on, 'baseline', 8); e_o = S.run(cfg, on, 'offset', 8)
            e_l = np.concatenate([one_trial_lsq(cfg, 1000 + t, on) for t in range(8)])
            say('%6d %12.2f %12.2f %12.2f' % (K, S.stats(e_b)['rms'], S.stats(e_o)['rms'], S.stats(e_l)['rms']))
        except Exception as ex:
            say('%6d layout failed: %s' % (K, ex))
    say('   (%.0f s)' % (time.time() - t0))
    os.makedirs('out', exist_ok=True)
    open('out/s58_decorrelator.txt', 'w').write('\n'.join(lines) + '\n')
    print('saved out/s58_decorrelator.txt')
