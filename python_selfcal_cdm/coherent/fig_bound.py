"""Fig. 10 of the paper: the worst-case bound (eq:bound) checked against the power model (a) and the coherent model (b).

Reads the capR_* records of campaign.py (all reflection orders, R = 1, 3, 10, 30 %, K = 4..32,
two random layouts each, 300-MHz line, one source realization) and, for the linewidth comparison,
capRb_* (1 GHz), capRbN_* (1 GHz, no detector noise) and capRbL_* (10 GHz). Records are decoded
without references or drift. The common chirp offset of each campaign (mean error of its R = 1 %,
K = 4 records, where every other mechanism is below 5 pm) is subtracted, as a reference would do.
For every grating the bound
    bound_k = 1.5 sum_j |Rule A(Delta_jk)| + 0.86 a_g sigma + 0.86 (K-1) sigma / (N T_k)
is evaluated from the stored layout: detunings give the Rule A sum and T_k (Gaussian lines of the
nominal width), positions give the ghost collisions a_g with a triangular overlap over one chip.

Findings (campaign.py capRd re-ran the R = 1 % and 3 % layouts without multiple reflections):
ghosts change the errors by at most 5 pm, detector noise changes nothing (capRbN), and the excess
over the bound falls with the linewidth, so it is the beat noise between the returns. RMS excess at
R = 1 %, K = 16: 16 pm (300 MHz), 4 pm (1 GHz), 1.4 pm (10 GHz).

Figure: (a) power model (../out/s55_results.npz from ../s55_error_bound.py),
        (b) coherent model at 300 MHz (open) and 1 GHz (filled), both as error of every grating against its bound.
Without the raw records the figure is drawn from cache/bound_field.npz (the decoded errors and bounds).
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ZP = ROOT
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ZP))
import decode as D
import figstyle as FS
FS.apply()
plt.rcParams.update({'pdf.fonttype': 42})

CACHE = HERE / 'records'
OUT = ROOT / 'figs'
SUMMARY = HERE / 'cache' / 'bound_field.npz'
SIG_PM = 250.0 / 2.35482
CMAX = np.sqrt(2.0) * np.exp(-0.5)
C0 = 299792458.0
N_CHIPS = 127


def rule_a_sum(det, R):
    K = len(det)
    out = np.zeros(K)
    for k in range(K):
        d = (det[:k] - det[k]) / SIG_PM
        out[k] = np.sum(np.abs((4.0 / 3.0) * np.sqrt(2.0 / 3.0) * R * d * np.exp(-d * d / 3.0))) * SIG_PM
    return out


def transmission_to(det, R):
    """Two-pass transmission of the upstream gratings at the Bragg wavelength of every grating."""
    K = len(det)
    T = np.ones(K)
    for k in range(K):
        for j in range(k):
            T[k] *= (1.0 - R * np.exp(-0.5 * ((det[j] - det[k]) / SIG_PM) ** 2)) ** 2
    return T


def ghost_ratio(z, det, R, dz):
    """Sum over third-order ghosts landing within one chip of grating k of their amplitude
    relative to the direct return of k (Gaussian lines, triangular overlap in delay)."""
    K = len(z)
    T = transmission_to(det, R)
    ratio = np.zeros(K)
    for a in range(K):
        for b in range(K):
            for c in range(K):
                if not (b < a and b < c):
                    continue
                zg = z[a] - z[b] + z[c]
                cbar = (det[a] + det[b] + det[c]) / 3.0
                spread = (det[a] - cbar) ** 2 + (det[b] - cbar) ** 2 + (det[c] - cbar) ** 2
                amp = R ** 3 * np.exp(-0.5 * spread / SIG_PM ** 2) * T[min(a, b, c)]
                for k in range(K):
                    w = 1.0 - abs(zg - z[k]) / dz
                    if w > 0:
                        g = amp * np.exp(-0.5 * ((cbar - det[k]) / (SIG_PM / np.sqrt(3.0))) ** 2)
                        ratio[k] += w * g / (R * T[k])
    return ratio


def record(prefix, R, K, t):
    f = CACHE / ('%s_R%02d_K%03d_t%d.npz' % (prefix, int(R * 100), K, t))
    if not f.exists():
        return None
    z = np.load(f, allow_pickle=True)
    e = D.analyze(z, refs=False, drift=False, stab=0)[0]
    det = z['det_pm']
    dz = C0 / (2.0 * float(z['n_group']) * float(z['chip_rate']))
    T = transmission_to(det, R)
    b = (1.5 * rule_a_sum(det, R) + CMAX * ghost_ratio(z['z'], det, R, dz) * SIG_PM
         + CMAX * (K - 1) * SIG_PM / (N_CHIPS * T))
    return e, b


def check(prefix, Rs=(0.01, 0.03, 0.10), Ks=(4, 8, 16, 32), trials=(0, 1)):
    base = [record(prefix, 0.01, 4, t) for t in trials]
    bias = float(np.mean([x[0].mean() for x in base if x is not None]))
    err, bnd, RR, KK = [], [], [], []
    print('== %s: chirp offset %.1f pm' % (prefix, bias))
    for R in Rs:
        for K in Ks:
            for t in trials:
                r = record(prefix, R, K, t)
                if r is None:
                    continue
                err.append(np.abs(r[0] - bias)); bnd.append(r[1]); RR.append(np.full(K, R)); KK.append(np.full(K, K))
    err, bnd, RR, KK = map(np.concatenate, (err, bnd, RR, KK))
    for R in Rs:
        for K in Ks:
            m = np.isclose(RR, R) & (KK == K)
            if m.any():
                exc = np.maximum(0.0, err[m] - bnd[m])
                print('   R=%3.0f%% K=%2d  sd %5.1f  max|e| %6.1f  max bound %9.1f  violations %5.1f%%  rms excess %5.1f'
                      % (100 * R, K, err[m].std(), err[m].max(), bnd[m].max(), 100 * np.mean(err[m] > bnd[m] + 0.05), np.sqrt(np.mean(exc ** 2))))
    return err, bnd, RR, KK


if __name__ == '__main__':
    res = {p: check(p) for p in ('capR', 'capRb', 'capRbN', 'capRbL') if any(CACHE.glob(p + '_R*_K*_t*.npz'))}   # only the sets that exist in the cache
    from_records = bool(res)
    if not from_records:                    # no raw records: decoded errors and bounds from the cache
        zc = np.load(SUMMARY)
        res = {p: tuple(zc['%s_%s' % (p, n)] for n in ('err', 'bound', 'R', 'K')) for p in ('capR', 'capRb', 'capRbN', 'capRbL') if ('%s_err' % p) in zc.files}
        print('records absent, figure drawn from', SUMMARY.name)
    pw = np.load(ZP / 'out/s55_results.npz')

    cols = {0.01: FS.C_GOOD, 0.03: FS.GREEN, 0.10: FS.C_MEAS}
    mk = {0.01: 's', 0.03: '^', 0.10: 'o'}
    fig, ax = plt.subplots(1, 2, figsize=(3.5, 1.9), layout='constrained')
    lim = (0.1, 500)
    for a, (err, bnd, RR, KK) in ((ax[0], (pw['err'], pw['bound'], pw['R'], pw['K'])), (ax[1], res['capR'])):
        for R in (0.01, 0.03, 0.10):
            m = np.isclose(RR, R)
            a.loglog(bnd[m], err[m], mk[R], ms=2.2, mfc='none', mec=cols[R], mew=0.6, alpha=0.7,
                     label='$R_0=%g\\%%$' % (100 * R))
        if a is ax[1] and 'capRb' in res:      # 1-GHz line: filled markers of the same shape
            errb, bndb, RRb, KKb = res['capRb']
            for R in (0.01, 0.03, 0.10):
                m = np.isclose(RRb, R)
                a.loglog(bndb[m], errb[m], mk[R], ms=2.0, mfc=cols[R], mec=cols[R], mew=0.4, alpha=0.6)
        a.plot(lim, lim, '-', color=FS.C_THEORY, lw=0.9)
        a.set_xlim(lim); a.set_ylim(lim)
        a.set_aspect('equal', adjustable='box')
        a.set_xlabel('Bound $\\delta\\lambda_k^{\\max}$ [pm]')
        a.legend(fontsize=5.5, loc='upper left', handlelength=1.2, borderaxespad=0.3, labelspacing=0.2)
    ax[0].set_ylabel('Error [pm]')
    ax[1].text(0.97, 0.05, 'open: 300 MHz\nfilled: 1 GHz', transform=ax[1].transAxes, ha='right', va='bottom', fontsize=5.2, color='0.3')
    ax[1].set_yticklabels([])
    FS.letter(ax[0], 'a'); FS.letter(ax[1], 'b')
    for ext in ('png', 'pdf'):
        try:
            fig.savefig(OUT / ('fig_bound_check.' + ext), dpi=300)
        except PermissionError:
            print('LOCKED, not written:', OUT / ('fig_bound_check.' + ext))
    if from_records:
        np.savez(SUMMARY, **{'%s_%s' % (p, n): v for p, r in res.items() for n, v in zip(('err', 'bound', 'R', 'K'), r)})
    print('saved', OUT / 'fig_bound_check.pdf')
