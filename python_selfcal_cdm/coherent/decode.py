"""decode.py - readout of the records of s70_field.run_array with the chain of the paper:
12-bit quantization, mean removed, circular correlation with the bipolar replica, read at the nearest
sample of the nominal delay, Gaussian fit of the peak (common.gauss_fit_peak), optional sequential
deshadowing (peel), reference correction with a drifted tuning table and a 1-pm reference stability.
"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import common as C

C0 = 299792458.0
ADC_BITS = 12
DRIFT_OFF_PM, DRIFT_GAIN, REF_STAB_PM = 20.0, 0.003, 1.0
DRIFT_V0 = 0.030                                            # V, offset voltage acting through dlambda/dV
# measured tuning curve of the BW10 VCSEL (quartic fit, nm vs V), as in s18_source of the paper
TUNE_P4 = np.array([-7.14581756e-05, 1.94718096e-03, -4.51295135e-02, -9.75634596e-02, 1.56974137e+03])
TUNE_V = np.linspace(0.0, 14.0, 2801)
TUNE_L = np.polyval(TUNE_P4, TUNE_V)
TUNE_C = 0.5 * (TUNE_L[0] + TUNE_L[-1])


def axis_error(x_pm, drift=1.0):
    """Reported minus true wavelength [pm] of a drifted tuning table at sweep position x_pm (pm from the
    sweep centre): 20-pm offset, 0.3 percent gain from V = 0, and 30 mV of offset voltage through the local
    slope of the measured curve. Same model as the chirp and drift section of the paper."""
    lam = TUNE_C + np.asarray(x_pm, float) / 1000.0
    V = np.interp(lam, TUNE_L[::-1], TUNE_V[::-1])
    lam_t = np.polyval(TUNE_P4, V)
    slope = np.polyval(np.polyder(TUNE_P4), V)
    return drift * (DRIFT_OFF_PM + DRIFT_GAIN * (lam_t - TUNE_L[0]) * 1000.0 + DRIFT_V0 * slope * 1000.0)


def correlate(y, code, spc, mode="baseline", adc=True):
    """Delay profile: quantization, mean removed, circular correlation with the replica 2c-1."""
    n = len(code) * spc
    y = np.asarray(y, float)[:n]
    if adc:
        fsr = 1.2 * max(y.max(), 1e-12)
        y = np.round(y / fsr * 2 ** (ADC_BITS - 1)) / 2 ** (ADC_BITS - 1) * fsr
    y = y - y.mean()
    c = np.repeat(2.0 * np.asarray(code) - 1.0, spc)
    X = np.fft.ifft(np.fft.fft(y) * np.conj(np.fft.fft(c))).real / (len(code) * spc)
    if mode == "offset":
        X = X - np.median(X)
    return X


def samples_of(z, ng, chip_rate, spc):
    """Delay of the return from distance z [m] in samples."""
    return 2.0 * z * ng / C0 * chip_rate * spc


def find_offset(X, bin_nom, spc):
    """Shift of the peak from the nominal bin (delay of the receiver filter), in samples."""
    n = X.shape[-1]
    win = (np.arange(int(round(bin_nom)) - spc // 2, int(round(bin_nom)) + spc // 2)) % n
    prof = X[:, win].max(axis=0) if X.ndim == 2 else X[win]
    return int(win[np.argmax(prof)]) - int(round(bin_nom))


def read_whole_chip(X, bin_samples, off, spc):
    """Read at the nearest WHOLE chip (an earlier decoder). At 100 Mchip/s and 4-m spacings (3.92 chips)
    this leaves the top of the correlation peak, hence read_nearest below."""
    n = X.shape[-1]
    idx = (int(round(bin_samples / spc)) * spc + off) % n
    return X[:, idx]


def read_nearest(X, bin_samples, off, spc):
    """Read at the nearest SAMPLE of the nominal delay (calibrated delay plus filter delay)."""
    n = X.shape[-1]
    idx = (int(round(bin_samples)) + off) % n
    return X[:, idx]


def read_exact(X, bin_samples, off):
    n = X.shape[-1]
    x = (bin_samples + off) % n
    i0 = int(np.floor(x)); f = x - i0
    return (1 - f) * X[:, i0 % n] + f * X[:, (i0 + 1) % n]


def peel(S, x, floor=0.05, R0=None):
    """Sequential deshadowing. R0: known reflectivities of the gratings (design values), used for the
    transmission instead of the fitted amplitude when given."""
    K = S.shape[0]
    T = np.ones(len(x)); cent = np.empty(K)
    for k in range(K):
        y = S[k] / np.maximum(T, floor)
        amp, mu, sg, _ = C.gauss_fit_full(x, y)
        cent[k] = mu
        a = float(R0[k]) if R0 is not None else np.clip(amp, 0, 0.99)
        T = T * (1.0 - a * np.exp(-0.5 * ((x - mu) / max(sg, 1e-3)) ** 2)) ** 2
    return cent


def analyze(z, mode="baseline", do_peel=False, refs=True, drift=True, stab=REF_STAB_PM, nearest=True,
            adc=True, seed_offset=0, ref_nub_pm=None, verbose=False, fit_order=1, peel_R0=True):
    """z: npz written by s70_field.run_array. Returns the errors [seeds x K] in pm (against the nominal detuning)."""
    code = z["code"]; spc = int(z["spc"]); rate = float(z["chip_rate"]); ng = float(z["n_group"])
    x = z["x_pm"]; det = z["det"]
    zs, ds, zr, dr = z["z"], z["det_pm"], z["zr"], z["detr_pm"]
    K, J = len(zs), len(zr)
    bins_s = samples_of(zs, ng, rate, spc); bins_r = samples_of(zr, ng, rate, spc)
    errs = []
    for si in range(det.shape[0]):
        X = np.array([correlate(det[si, m], code, spc, mode, adc) for m in range(len(x))])
        # filter delay from the first element (a reference when present, otherwise sensor 1)
        b0 = bins_r[0] if J else bins_s[0]
        off = find_offset(X, b0, spc)
        rd = (lambda b: read_nearest(X, b, off, spc)) if nearest else (lambda b: read_exact(X, b, off))
        S = np.array([rd(b) for b in bins_s]); Sr = np.array([rd(b) for b in bins_r])
        # scale: peel needs amplitudes in units of R, normalized to the mean reference amplitude
        # (references have a known R and are isolated), otherwise to the maximum
        scale = np.mean([Sr[j].max() / float(z["Rr"][j]) for j in range(J)]) if J else S.max() / float(z["R"][0])
        S = S / scale; Sr = Sr / scale
        if do_peel:
            cent = peel(S, x, R0=(z["R"] if peel_R0 else None))
        else:
            cent = np.array([C.gauss_fit_peak(x, S[k]) for k in range(K)])
        rng = np.random.default_rng(1000 + si + seed_offset)
        dr_ax = axis_error(x) if drift else np.zeros_like(x)
        cent_rep = cent - np.interp(cent, x, dr_ax)
        if J and refs:
            # references have known wavelengths: fit within +-250 pm of the nominal value (immune to beat noise)
            cent_r = []
            for j in range(J):
                c = C.gauss_fit_peak(x, Sr[j])
                if abs(c - dr[j]) > 40.0:          # wrong peak (beat noise): fit within a window around the nominal value
                    w = np.abs(x - dr[j]) <= 250.0
                    c = C.gauss_fit_peak(x[w], Sr[j][w])
                cent_r.append(c)
            cent_r = np.array(cent_r)
            cent_r_rep = cent_r - np.interp(cent_r, x, dr_ax)
            nub_r = np.asarray(ref_nub_pm if ref_nub_pm is not None else dr, float)
            err_r = cent_r_rep - (nub_r + rng.normal(0.0, stab, J))
            p = np.polyfit(nub_r, err_r, min(fit_order, J - 1))    # least squares, order of the drift model
            cent_rep = cent_rep - np.polyval(p, cent_rep)
        errs.append(cent_rep - ds)
        if verbose:
            print("   seed %d: off %d, centres %s" % (si, off, np.array2string(cent, precision=1)))
    return np.array(errs)


def stats(e):
    a = np.abs(e)
    return dict(rms=float(np.sqrt(np.mean(e ** 2))), p90=float(np.percentile(a, 90)), mx=float(a.max()),
                per=np.sqrt(np.mean(e ** 2, axis=0)))
