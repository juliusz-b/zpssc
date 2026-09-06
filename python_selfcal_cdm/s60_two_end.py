"""s60_two_end.py - capacity of one band when the array is read from both
ends and the receiver decorrelates: is K set by the delay bins alone?

The two ideas of s58 and s59 together, on the full chain of s56 (optics with
shadowing and every third-order ghost, chirp kernel, receiver noise, filter,
ADC). Four readouts of the same array, same random draws:

  MF          matched filter from the source end, no correction (baseline)
  MF + seq    matched filter, sequential deshadowing with fitted transmissions
  DC + seq    decorrelating receiver (s58), sequential deshadowing
  DC two-end  decorrelating receiver from both ends, closed-form recovery
              R_k/(1-R_k) = sqrt(S_k^f S_k^b) / T_tot (s59), T_tot measured
              in transmission with the same relative noise

Strong gratings on purpose: R = 10 percent, 250-pm lines, one band of
+-200 pm, N = 127 at 100 Mchip/s (one chip = 1.02 m), gratings in distinct
chip bins with a random fraction of a chip, K from 8 to 96. References on a
stub, drift and 1-pm stability as in s56, chirp span 0.2 FWHM. The backward
reading uses the mirrored layout, so its ghosts and shadowing are its own.

Output: out/s60_two_end.txt
"""
import os, time, warnings
import numpy as np
import common as C
import s56_budget_e2e as S
import s58_decorrelator as D
warnings.filterwarnings('ignore')

NT = int(os.environ.get('S60_NTRIALS', '8'))


def fit_all(Sk):
    return np.array([C.gauss_fit_peak(S.LAM, Sk[k]) for k in range(Sk.shape[0])])


def correct_refs(cent, cent_r, nub_r, nref, drift, rng, on):
    cent_rep = cent - np.interp(cent, S.LAM, drift)
    cent_r_rep = cent_r - np.interp(cent_r, S.LAM, drift)
    if on['refs']:
        stab = rng.normal(0.0, S.REF_STAB_PM / S.PM, nref)
        err_r = cent_r_rep - (nub_r + stab)
        p = np.polyfit(nub_r, err_r, nref - 1) if nref > 1 else np.array([err_r.mean()])
        cent_rep = cent_rep - np.polyval(p, cent_rep)
    return cent_rep


def one_trial(cfg, seed, on):
    rng = np.random.default_rng(seed)
    K = cfg['K']
    tau_s, tau_r = S.make_layout(cfg, rng)
    nub = rng.uniform(-S.W_GHZ, S.W_GHZ, K)
    asym = rng.uniform(-S.ASYM_SPREAD, S.ASYM_SPREAD, K)
    nub_r = np.linspace(-0.9 * S.W_GHZ, 0.9 * S.W_GHZ, cfg['nref'])
    F = cfg['fwhm'] / S.PM
    std_chirp = cfg['chirp'] * F / 2.0
    off_fn = lambda nb: std_chirp * (0.20 + 0.20 * nb / S.W_GHZ + 0.15 * (nb / S.W_GHZ) ** 2 + 0.10 * np.sin(2.5 * nb / S.W_GHZ))
    drift = (S.DRIFT_OFF_PM + S.DRIFT_GAIN * S.LAM * S.PM) / S.PM if on['drift'] else np.zeros(S.M_STEPS)
    known = np.concatenate([tau_s, tau_r])
    # forward
    delays, A = S.spectra_matrix(cfg, nub, asym, tau_s, tau_r, nub_r, on, std_chirp, off_fn)
    prof = S.receiver(cfg, delays, A, on, np.random.default_rng(seed + 7), 'baseline')
    S_mf = np.array([S.read(prof, t, False) for t in tau_s]); Sr_mf = np.array([S.read(prof, t, False) for t in tau_r])
    est = D.receiver_lsq(cfg, delays, A, on, np.random.default_rng(seed + 7), known)
    S_f = est[:, :K].T; Sr_f = est[:, K:].T
    # backward: mirrored positions, gratings renumbered from the far end
    tau_b = tau_s.max() + tau_s.min() - tau_s          # same spacings, reversed order
    order = np.argsort(tau_b)
    delays_b, A_b = S.spectra_matrix(cfg, nub[order], asym[order], tau_b[order], tau_r, nub_r, on, std_chirp, off_fn)
    est_b = D.receiver_lsq(cfg, delays_b, A_b, on, np.random.default_rng(seed + 11), np.concatenate([tau_b[order], tau_r]))
    S_b = np.empty_like(S_f); S_b[order] = est_b[:, :K].T
    # transmission of the whole array, measured once with the same relative noise
    R = cfg['R']
    Rk = np.array([np.clip(R * np.exp(-0.5 * ((S.LAM - nub[k]) / (F / 2.35482)) ** 2) * (1.0 + asym[k] * np.tanh((S.LAM - nub[k]) / (F / 2.35482))), 0, None) for k in range(K)])
    Ttot = np.prod(1.0 - Rk, axis=0)
    if on['noise']:
        Ttot = Ttot * (1.0 + rng.normal(0.0, 3e-4, Ttot.shape))
    q = np.sqrt(np.clip(S_f * S_b, 0, None)) / np.clip(Ttot, 1e-6, None)[None, :]
    R_bi = q / (1.0 + q)
    # centres
    truth = np.array([C.gauss_fit_peak(S.LAM, Rk[k]) for k in range(K)])
    cr_mf = fit_all(Sr_mf); cr_f = fit_all(Sr_f)
    out = {}
    out['MF'] = correct_refs(fit_all(S_mf), cr_mf, nub_r, cfg['nref'], drift, np.random.default_rng(seed + 1), on)
    out['MF + seq'] = correct_refs(S.peel(S_mf, nub), cr_mf, nub_r, cfg['nref'], drift, np.random.default_rng(seed + 1), on)
    out['DC + seq'] = correct_refs(S.peel(S_f, nub), cr_f, nub_r, cfg['nref'], drift, np.random.default_rng(seed + 1), on)
    out['DC two-end'] = correct_refs(fit_all(R_bi), cr_f, nub_r, cfg['nref'], drift, np.random.default_rng(seed + 1), on)
    return {k: (v - truth) * S.PM for k, v in out.items()}


if __name__ == '__main__':
    lines = []
    say = lambda s='': (print(s), lines.append(s))
    on = {k: True for k in S.ALL}
    modes = ['MF', 'MF + seq', 'DC + seq', 'DC two-end']
    t0 = time.time()
    say('R = 10 percent, N = 127 at 100 Mchip/s, distinct chip bins, chirp 0.2 FWHM, %d trials. RMS Bragg error [pm], p90 in brackets' % NT)
    say('%4s ' % 'K' + ''.join('%22s' % m for m in modes))
    for K in (8, 16, 32, 48, 64, 80, 96):
        cfg = dict(K=K, R=0.10, fwhm=250.0, nbits=7, B=100e6, spacing='bins', d_min=1.0, nref=3, peel=True, chirp=0.2)
        res = {m: [] for m in modes}
        for t in range(NT):
            r = one_trial(cfg, 3000 + t, on)
            for m in modes:
                res[m].append(r[m])
        row = '%4d ' % K
        for m in modes:
            e = np.concatenate(res[m]); row += '%13.2f (%6.2f)' % (np.sqrt(np.mean(e ** 2)), np.percentile(np.abs(e), 90))
        say(row + '   %.0f s' % (time.time() - t0))
    os.makedirs('out', exist_ok=True)
    open('out/s60_two_end.txt', 'w').write('\n'.join(lines) + '\n')
    print('saved out/s60_two_end.txt')
