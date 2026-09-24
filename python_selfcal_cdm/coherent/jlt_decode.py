"""Decoding of the 50-sensor records (records/jlt/jltc_*.npz, written by jlt_chain.py):
nearest-sample read of every delay bin, then (a) the Gaussian fit over the +-S = 450 pm window around the NOMINAL
wavelength (all samples of the window, free baseline, width bounded to 60..600 pm, the estimator used in the paper:
err, errc) and (b) the centroid of [JLT] eq. (22) restricted to the same window (err_c).
The local fit of Table III (samples above 30 % of the peak) is not used here, because the 80-pm steps of [JLT] leave
only three samples above half maximum and that fit fails on the flat returns of the co-tuned array (sigma > 100 pm). The nominal wavelength is the design value, not the true one: for the detuned
case D4x it is the band layout (det_bands of jlt_chain.py). Only for case 2 of [JLT] (J5mb, shifts up to 0.58 nm,
beyond the window) the window follows the sensors, as a tracking receiver would. reference correction of offset and slope where the record has references
(the same estimator applied to the references), otherwise the common offset removed.
Writes cache/jlt/jltc_<name>_dec.npz (S, det, z, err, errc, err_c) for fig_jlt.py and prints a table."""
import sys, os, glob, numpy as np
from scipy.optimize import curve_fit
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import decode as D
import common as C
T = os.path.join(HERE, 'records', 'jlt')
TC = os.path.join(HERE, 'cache', 'jlt'); os.makedirs(TC, exist_ok=True)
S_PM = 450.0          # operating range +-200 pm plus the 250-pm line, at 80-pm steps
K_JLT = 50
_band = np.repeat([0, 1, 2], [17, 17, 16]); DET_BANDS = np.zeros(K_JLT)
for _b in range(3):
    _idx = np.flatnonzero(_band == _b); DET_BANDS[_idx] = -800.0 + 800.0 * _b + np.linspace(-350.0, 350.0, len(_idx))
CHIP = {'D0': 511 / 16.352e-6, 'J5': 511 / 16.352e-6, 'J50': 511 / 16.352e-6, 'J5m': 511 / 16.352e-6, 'J5m50': 50e6, 'J5mb': 511 / 16.352e-6, 'D4x': 50e6, 'D4y': 50e6, 'D4g': 50e6, 'B2': 50e6, 'B5': 50e6, 'D4j': 50e6, 'D1': 511 / 16.352e-6, 'D2': 50e6, 'D3': 50e6, 'D4': 50e6}


def near(x, y, x0, W=8):
    """Samples around the largest value within +-S of x0, the rest flattened."""
    m = np.abs(x - x0) <= S_PM
    y2 = np.where(m, y, y[m].min())
    i0 = int(np.argmax(y2))
    lo, hi = max(0, i0 - W), min(len(x), i0 + W + 1)
    return x[lo:hi], y2[lo:hi]


def _gauss(x, a, x0, s, b):
    return a * np.exp(-0.5 * ((x - x0) / s) ** 2) + b


def gauss_win(x, y, x0):
    """Gaussian fit to all samples within +-S_PM of the nominal wavelength, free baseline, width 60..600 pm.
    Returns nan when the fit does not converge (counted as an error of S by the caller)."""
    m = np.abs(x - x0) <= S_PM
    xs, ys = x[m], y[m]
    p0 = [ys.max(), xs[int(np.argmax(ys))], 150.0, 0.0]
    try:
        p, _ = curve_fit(_gauss, xs, ys, p0=p0, bounds=([0.0, x0 - S_PM, 60.0, -np.inf], [np.inf, x0 + S_PM, 600.0, np.inf]), maxfev=20000)
        return float(p[1])
    except Exception:
        return np.nan


def cen_win(x, y, x0):
    """Centroid (22) of [JLT] over the samples within +-S_PM of the nominal wavelength, negative values clipped."""
    m = np.abs(x - x0) <= S_PM
    xw, yw = x[m], np.clip(y[m], 0, None)
    return (xw * yw).sum() / yw.sum()


def decode(path):
    z = np.load(path, allow_pickle=True)
    name = os.path.basename(path)[5:-4]
    import re
    base = re.sub(r'L\d+$', '', re.sub(r'R\d+', '', re.sub(r'A(hn|hm|tk|ga|bl)', '', name.split('_')[0])).rstrip('r')).rstrip('w')
    base = base if base in CHIP else base[:2]
    chip = CHIP[base]
    code = z['code']; spc = int(z['spc']); ng = float(z['n_group']); x = z['x_pm']
    det = z['det_pm']; K = len(det)
    nom = DET_BANDS if base in ('D4x', 'D4y', 'D4g') else det          # window center: design wavelength (D4x), else the true one equals it (J5mb: tracking)
    bins = D.samples_of(z['z'], ng, chip, spc)
    X = np.array([D.correlate(z['det'][0, m], code, spc) for m in range(len(x))])
    n = X.shape[1]
    off = D.find_offset(X, bins[0], spc)
    S = np.array([X[:, (int(round(b)) + off) % n] for b in bins]); scale = S[0].max(); S = S / scale
    err_c = np.array([cen_win(x, S[k], nom[k]) for k in range(K)]) - det       # centroid within +-S_PM of the nominal wavelength, as (22) of [JLT]
    fit = np.array([gauss_win(x, S[k], nom[k]) for k in range(K)])
    fit = np.where(np.isfinite(fit), fit, nom + S_PM)      # a failed fit counts as an error of S
    err = fit - det
    if len(z['zr']):
        bins_r = D.samples_of(z['zr'], ng, chip, spc)
        Sr = np.array([X[:, (int(round(b)) + off) % n] for b in bins_r]) / scale
        fit_r = np.array([gauss_win(x, Sr[j], z['detr_pm'][j]) for j in range(len(bins_r))])
        fit_r = np.where(np.isfinite(fit_r), fit_r, z['detr_pm'] + S_PM)
        p = np.polyfit(z['detr_pm'], fit_r - z['detr_pm'], 1)
        cen_r = np.array([cen_win(x, Sr[j], z['detr_pm'][j]) for j in range(len(bins_r))])
        pc = np.polyfit(z['detr_pm'], cen_r - z['detr_pm'], 1)
        errc = err - np.polyval(p, det); err_cc = err_c - np.polyval(pc, det); refs = 2
    else:
        errc = err - err.mean(); err_cc = err_c - err_c.mean(); refs = 0
    np.savez(os.path.join(TC, 'jltc_%s_dec.npz' % name), x_pm=x, S=S, det=det, z=z['z'], err=err, errc=errc, err_c=err_cc)
    return name, refs, err_cc, errc, S.max(axis=1)


if __name__ == '__main__':
    rows = []
    for f in sorted(glob.glob(os.path.join(T, 'jltc_*.npz'))):
        if f.endswith('_dec.npz') or 'ref1200' in f:
            continue
        name, refs, ec, ef, amp = decode(f)
        print('%-8s refs %d | Gaussian fit sigma %6.1f RMS %6.1f max %6.1f | centroid sigma %6.1f RMS %6.1f max %6.1f pm | amp 10..50 %s' % (
            name, refs, ef.std(), np.sqrt(np.mean(ef ** 2)), np.abs(ef).max(), ec.std(), np.sqrt(np.mean(ec ** 2)), np.abs(ec).max(),
            np.round(amp[[9, 19, 29, 39, 49]], 2)), flush=True)
