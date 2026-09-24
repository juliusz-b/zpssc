"""Beat-noise floor of the correlogram against a closed-form estimate.

Estimate: for K returns of powers P_k sharing one wavelength step, unipolar code of N chips at chip rate B, source line
Delta_nu (Lorentzian), the cross terms of the detected intensity carry the power (1/4) sum_{k<l} 2 P_k P_l, of which the
chip integration passes the fraction f = (2/pi) arctan(B / (2 Delta_nu)). After correlation over N chips the RMS floor
relative to the peak of an isolated grating of power P_1 is

    floor = (1+eps)/(1-eps) sqrt( f [ (sum P)^2 - sum P^2 ] / N ) / P_1 ,

where eps is the power of a 0 chip relative to a 1 chip (finite extinction of the directly modulated laser, 0.30 in the
model at 33 +- 18 mA), which keeps the cross terms alive during the 0 chips and lowers the correlation peak.

Checked on one-step records (common Bragg wavelength, direct paths only, no detector noise, receiver 0.75 B) for the
[JLT] array (50 gratings 2.5 m apart, R = 1 %, P_k = P_1 (1-R)^(2(k-1))) at Delta_nu = 0.305, 1, 3 GHz and
B = 31.25, 100 Mchip/s, and for one band of 17 gratings. Field oversampling raised with the linewidth.
python beat_floor.py -> table on stdout, figs/fig_beat_floor_check.png, cache/beat_floor_rows.npy (Fig. 11b, fig_validation.py)"""
import sys, os, time, functools
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
import s70_field as F
import decode as D
import figstyle as FS
FS.apply()
F.photocurrent = functools.partial(F.photocurrent, resp=0.8, nep=15e-12)
T = str(HERE / 'records' / 'beat'); os.makedirs(T, exist_ok=True)
C0 = 299792458.0
R = 0.01
EPS = 0.30                     # 0-chip power relative to the 1-chip power (measured on the reference record)
mseq = F._mls(9); N = 511
CASES = []
for K in (50, 17):
    zs = np.sort(2.5 * np.arange(1, K + 1) + np.random.default_rng(11).uniform(-0.05, 0.05, K))
    for lw in (305e6, 1e9, 3e9):
        for B in ((511 / 16.352e-6, 100e6) if K == 50 else (511 / 16.352e-6,)):
            CASES.append((K, zs, lw, B))


def rec(name, zs, lw, B):
    out = '%s/beat_%s.npz' % (T, name)
    if not os.path.exists(out):
        osr = 4 if lw <= 1e9 else 8
        t0 = time.time()
        F.run_array([dict(z=float(z), det=0.0, R=R, g='bl250') for z in zs], [], code=mseq, chip_rate=B, x_pm=np.array([0.0]),
                    seeds=(1, 2), spc=64, linewidth_hz=lw, det_noise=False, ghosts=False, laser_params=dict(eta=0.063),
                    out=out, verbose=False, osr=osr)
        print(name, '%.0f s' % (time.time() - t0), flush=True)
    return np.load(out, allow_pickle=True)


def main():
    rows = []
    for K, zs, lw, B in CASES:
        name = 'K%d_lw%d_B%d' % (K, lw / 1e6, B / 1e6)
        ref = rec('ref_' + name, zs[:1], lw, B)
        z = rec(name, zs, lw, B)
        code = z['code']; spc = int(z['spc']); ng = float(z['n_group'])
        Lc = C0 / (2 * ng * B); zax = np.arange(len(code) * spc) / spc * Lc
        unit = np.mean([D.correlate(ref['det'][s, 0], code, spc).max() for s in range(2)])
        floors = []
        for s in range(2):
            X = D.correlate(z['det'][s, 0], code, spc) / unit
            empty = (zax > zs[-1] + 40) & (zax < 0.9 * zax.max())
            floors.append(X[empty].std())
        P = (1 - R) ** (2 * np.arange(K))
        f = (2 / np.pi) * np.arctan(B / (2 * lw))
        est = (1 + EPS) / (1 - EPS) * np.sqrt(f * (P.sum() ** 2 - (P ** 2).sum()) / N)
        rows.append((K, lw, B, np.mean(floors), est))
        print('K=%2d  line %5.0f MHz  B %6.2f Mchip/s | model floor %.3f | estimate %.3f | ratio %.2f' % (K, lw / 1e6, B / 1e6, np.mean(floors), est, np.mean(floors) / est), flush=True)
    rows = np.array(rows)
    fig, ax = plt.subplots(figsize=(3.5, 2.6), layout='constrained')
    for K, mk, col in ((50, 'o', FS.VERM), (17, 's', FS.BLUE)):
        m = rows[:, 0] == K
        ax.plot(rows[m, 4], rows[m, 3], mk, color=col, ms=5, mfc='none', label='K = %d' % K)
    lim = [0.02, 0.6]
    ax.plot(lim, lim, 'k--', lw=0.8, label='estimate')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel('Estimated floor'); ax.set_ylabel('Floor in the coherent model'); ax.legend(fontsize=6)
    png = HERE.parent / 'figs' / 'fig_beat_floor_check.png'
    fig.savefig(png, dpi=150); print('saved', png)
    np.save(HERE / 'cache' / 'beat_floor_rows.npy', rows)


if __name__ == '__main__':
    main()
