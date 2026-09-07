"""s70_field.py - field-level model of the coded FBG interrogator, fast and self-contained.

Gratings and fiber: coupled-mode transfer matrices (Erdogan 1997), cascaded with the fiber
propagation phase between gratings. The total reflection coefficient r(w) of the serial array
therefore contains shadowing, every order of multiple reflections, grating dispersion and the
interference between paths of equal optical length. References sit on a separate branch whose
field is added to the sensor return.

Source: single-mode rate equations (carrier density, photon density, phase with the linewidth
enhancement factor) integrated once per noise realization at a fine step and resampled to the
record grid. The intensity code, the adiabatic and transient chirp and the relaxation
oscillations follow from the equations. A Lorentzian linewidth is added as a random walk of the
phase, with a new realization at every sweep step, and an optional thermal chirp as a first-order
low-pass of the drive current.

Readout: the record is periodic, so the array acts as a linear filter: E_out = IFFT(r(w) FFT(E_in)).
The photodiode gives |E_out|^2 with shot and thermal noise, then a causal fourth-order Bessel
low-pass at 0.75 B, sampled at spc samples per chip. The output has the same layout as the
records of the VPIphotonics runs (det[seed, step, sample]), so the same decoding applies.

Typical cost: one sweep of 64 steps for 32 gratings in a few seconds.
"""
import os
import numpy as np
from scipy.signal import bessel, lfilter

C0 = 299792458.0
N_EFF = 1.45
N_GROUP = 1.468
LAM0 = 1552.524e-9
Q_E = 1.602e-19
H_NU = 6.626e-34 * C0 / LAM0


# ----------------------------------------------------------------------------- gratings
def blackman(x):
    """Apodization profile on x in [-0.5, 0.5], peak 1 at the center."""
    return (1.0 + 1.19 * np.cos(2 * np.pi * x) + 0.19 * np.cos(4 * np.pi * x)) / 2.38


def grating_matrix(lam, lam_b, length, kappa0, apod="blackman", sections=40):
    """Transfer matrix [R_in; S_in] = F [R_out; S_out] of one grating for every wavelength in lam.
    Returns F as arrays (f11, f12, f21, f22), each of lam.shape."""
    beta = 2 * np.pi * N_EFF / lam
    period = lam_b / (2 * N_EFF)
    sigma = beta - np.pi / period                       # detuning from the Bragg condition
    dz = length / sections
    x = (np.arange(sections) + 0.5) / sections - 0.5
    w = blackman(x) if apod == "blackman" else np.ones(sections)
    f11 = np.ones_like(lam, dtype=complex); f12 = np.zeros_like(f11); f21 = np.zeros_like(f11); f22 = np.ones_like(f11)
    for m in range(sections):
        k = kappa0 * w[m]
        g = np.sqrt((k ** 2 - sigma ** 2).astype(complex))
        g = np.where(np.abs(g) < 1e-12, 1e-12, g)
        ch, sh = np.cosh(g * dz), np.sinh(g * dz)
        a11 = ch - 1j * (sigma / g) * sh
        a12 = -1j * (k / g) * sh
        a21 = 1j * (k / g) * sh
        a22 = ch + 1j * (sigma / g) * sh
        n11 = f11 * a11 + f12 * a21; n12 = f11 * a12 + f12 * a22
        n21 = f21 * a11 + f22 * a21; n22 = f21 * a12 + f22 * a22
        f11, f12, f21, f22 = n11, n12, n21, n22
    return f11, f12, f21, f22


def grating_r(lam, lam_b, length, kappa0, **kw):
    f11, f12, f21, f22 = grating_matrix(lam, lam_b, length, kappa0, **kw)
    return f21 / f11


# calibrated Blackman gratings: (length, kappa0) giving the FWHM at the stated peak reflectivity
GRATINGS = {
    "bl250": dict(length=4.6e-3, fwhm=250.0),
    "bl100": dict(length=11.5e-3, fwhm=100.0),
    "uni40": dict(length=10.0e-3, fwhm=40.0, apod="uniform"),
}


def calibrate(gtype, R, lam_b=LAM0):
    """Find (length, kappa0) for the grating type so that peak reflectivity is R and FWHM matches.
    Cached per (gtype, R)."""
    key = (gtype, round(R, 5))
    if key in _CAL:
        return _CAL[key]
    g = GRATINGS[gtype]
    apod = g.get("apod", "blackman")
    lam = lam_b + np.linspace(-800e-12, 800e-12, 1601)

    def peak_and_fwhm(length, kappa0):
        Rr = np.abs(grating_r(lam, lam_b, length, kappa0, apod=apod)) ** 2
        half = Rr >= 0.5 * Rr.max()
        return Rr.max(), (lam[half].max() - lam[half].min()) * 1e12

    length = g["length"]
    for _ in range(3):
        lo, hi = 1.0, 3000.0
        for _ in range(40):                          # bisection on kappa0 for the peak
            mid = 0.5 * (lo + hi)
            if peak_and_fwhm(length, mid)[0] < R:
                lo = mid
            else:
                hi = mid
        kappa0 = 0.5 * (lo + hi)
        pk, fw = peak_and_fwhm(length, kappa0)
        if g.get("apod", "blackman") == "uniform":
            break
        length = length * fw / g["fwhm"]             # FWHM scales roughly as 1/length
    _CAL[key] = (length, kappa0, pk, fw)
    return _CAL[key]


_CAL = {}
_LCACHE = {}


# ----------------------------------------------------------------------------- array
def array_reflection_direct(lam, elements):
    """Field of the direct paths only (reflection of k through the transmissions of j < k, twice):
    the power-model picture without multiple reflections, for ablation against array_reflection."""
    order = np.argsort([e["z"] for e in elements])
    elems = [elements[i] for i in order]
    w = 2 * np.pi * C0 / lam
    w0 = 2 * np.pi * C0 / LAM0
    beta_f = (N_EFF * w0 + N_GROUP * (w - w0)) / C0
    r_tot = np.zeros_like(lam, dtype=complex)
    t_acc = np.ones_like(lam, dtype=complex)
    zprev = 0.0
    ph_acc = np.ones_like(lam, dtype=complex)
    for e in elems:
        ph_acc = ph_acc * np.exp(-1j * beta_f * (e["z"] - zprev))       # same sign as T21/T11 of the cascade
        f11, f12, f21, f22 = grating_matrix(lam, e["lam_b"], e["length"], e["kappa0"], apod=e.get("apod", "blackman"))
        r_k, t_k = f21 / f11, 1.0 / f11
        r_tot = r_tot + ph_acc ** 2 * t_acc ** 2 * r_k
        t_acc = t_acc * t_k
        zprev = e["z"]
    return r_tot


def array_reflection(lam, elements, lead_m=0.0):
    """Total complex reflection of a serial array. elements: list of dict(z, lam_b, length, kappa0, apod).
    z in meters from the circulator (sorted). Includes all orders of multiple reflections."""
    order = np.argsort([e["z"] for e in elements])
    elems = [elements[i] for i in order]
    w = 2 * np.pi * C0 / lam
    w0 = 2 * np.pi * C0 / LAM0
    beta_f = (N_EFF * w0 + N_GROUP * (w - w0)) / C0    # fiber propagation constant with the group delay
    T11 = np.ones_like(lam, dtype=complex); T12 = np.zeros_like(T11); T21 = np.zeros_like(T11); T22 = np.ones_like(T11)

    def mul(a11, a12, a21, a22):
        nonlocal T11, T12, T21, T22
        n11 = T11 * a11 + T12 * a21; n12 = T11 * a12 + T12 * a22
        n21 = T21 * a11 + T22 * a21; n22 = T21 * a12 + T22 * a22
        T11, T12, T21, T22 = n11, n12, n21, n22

    zprev = 0.0
    for e in elems:
        d = e["z"] - zprev + (lead_m if zprev == 0.0 else 0.0)
        ph = np.exp(1j * beta_f * d)
        mul(ph, 0.0, 0.0, 1.0 / ph)                         # fiber section
        mul(*grating_matrix(lam, e["lam_b"], e["length"], e["kappa0"], apod=e.get("apod", "blackman")))
        zprev = e["z"]
    return T21 / T11


# ----------------------------------------------------------------------------- laser
LASER = dict(
    tau_n=4.0e-9,          # carrier lifetime [s]
    tau_p=2.0e-12,         # photon lifetime [s]
    gamma=0.3,             # confinement factor
    vg=C0 / 4.0,           # group velocity in the cavity [m/s]
    a=3.3e-20,             # differential gain [m^2]
    n_tr=1.5e24,           # transparency carrier density [1/m^3]
    eps=3.0e-23,           # gain compression [m^3]
    beta_sp=1e-4,          # spontaneous emission factor
    alpha=4.0,             # linewidth enhancement factor
    volume=1.8e-16,        # active volume [m^3] (300 um x 3 um x 0.2 um), threshold about 15 mA
    eta=0.3,               # power per photon density scale (arbitrary units)
)


def laser_record(current, fs, seed=0, params=None, dt_int=5e-12, thermal=None):
    """Integrate the rate equations for the drive current [A] sampled at fs and return the complex
    field envelope sqrt(P) exp(i phi) on the same grid. thermal = (tau_s, hz_per_amp) adds a thermal
    frequency shift that follows the current through a first-order low-pass."""
    p = dict(LASER, **(params or {}))
    n = current.size
    sub = max(1, int(round(1.0 / (fs * dt_int))))
    dt = 1.0 / (fs * sub)
    I = np.repeat(current, sub)
    if thermal is not None:
        tau_s, hz_per_amp = thermal
        af = dt / tau_s
        Ith = np.empty_like(I); acc = I[0]
        for i in range(I.size):
            acc += af * (I[i] - acc); Ith[i] = acc
        dphi_th = 2 * np.pi * hz_per_amp * (Ith - Ith.mean())
    else:
        dphi_th = np.zeros_like(I)
    V, G, vg, a, ntr, eps, tn, tp, bsp, al = (p[k] for k in ("volume", "gamma", "vg", "a", "n_tr", "eps", "tau_n", "tau_p", "beta_sp", "alpha"))
    N, S, phi = ntr * 1.2, 1e18, 0.0
    n_th = ntr + 1.0 / (G * vg * a * tp)                 # threshold carrier density (uncompressed gain)
    # settle at the first current value
    for _ in range(int(20e-9 / dt)):
        g = vg * a * (N - ntr) / (1 + eps * S)
        dN = I[0] / (Q_E * V) - N / tn - g * S
        dS = G * g * S - S / tp + G * bsp * N / tn
        N += dt * dN; S = max(S + dt * dS, 1e10)
    # two periods of the drive, the second one is kept: periodic steady state at the record boundary
    out_s = np.empty(I.size); out_phi = np.empty(I.size)
    for rep in range(2):
        for i in range(I.size):
            g = vg * a * (N - ntr) / (1 + eps * S)
            dN = I[i] / (Q_E * V) - N / tn - g * S
            dS = G * g * S - S / tp + G * bsp * N / tn
            # index change follows the carrier density: chirp = alpha/2 * Gamma vg a (N - N_th),
            # which with gain compression gives the adiabatic term and with dN/dt the transient term
            dphi = 0.5 * al * G * vg * a * (N - n_th) + dphi_th[i]
            N += dt * dN; S = max(S + dt * dS, 1e10); phi += dt * dphi
            if rep == 1:
                out_s[i] = S; out_phi[i] = phi
    # mean emission frequency (power weighted) becomes the band center: removed here, added back in the readout
    dphi = np.diff(out_phi) / dt
    f_c = float(np.sum(out_s[1:] * dphi) / np.sum(out_s[1:]) / (2 * np.pi))
    t = np.arange(I.size) * dt
    out = np.sqrt(out_s) * np.exp(1j * (out_phi - 2 * np.pi * f_c * t))
    E = out.reshape(n, sub).mean(axis=1)
    return E * np.sqrt(p["eta"] / (np.mean(np.abs(E) ** 2) + 1e-30)), f_c      # mean |E|^2 = eta, arbitrary power unit


def drive_current(code, spc, bias, amp, rise_frac=0.25):
    """NRZ current with linear edges of rise_frac chip."""
    c = np.repeat(np.asarray(code, float), spc)
    k = max(1, int(round(rise_frac * spc)))
    ker = np.ones(k) / k
    c = np.convolve(np.r_[c[-k:], c, c[:k]], ker, mode="same")[k:-k]
    return bias + amp * (c - 0.5)


def phase_walk(n, fs, linewidth_hz, rng):
    if not linewidth_hz:
        return np.zeros(n)
    dphi = rng.normal(0.0, np.sqrt(2 * np.pi * linewidth_hz / fs), n)
    return np.cumsum(dphi)


# ----------------------------------------------------------------------------- readout
def photocurrent(E_out, fs, resp=1.0, p_scale=1e-3, nep=0.5e-12, det_noise=True, rng=None):
    P = p_scale * np.abs(E_out) ** 2
    i = resp * P
    if det_noise:
        rng = rng or np.random.default_rng()
        i = i + rng.normal(0.0, np.sqrt(2 * Q_E * resp * max(P.mean(), 0) * fs / 2), i.size)
        i = i + rng.normal(0.0, resp * nep * np.sqrt(fs / 2), i.size)
    return i


def run_array(sensors, refs=(), code=None, chip_rate=25e6, x_pm=None, seeds=(1,), spc=64, bias=0.033, amp=0.018,
              rise_frac=0.25, linewidth_hz=30e6, thermal=None, det_noise=True, laser_params=None, out=None, verbose=False, ghosts=True):
    """Sweep the laser over x_pm (pm around LAM0) and return records det[seed, step, sample] plus metadata,
    in the layout of the VPI runs. sensors, refs: dict(z, det, R, g)."""
    code = np.asarray(code if code is not None else _mls(7), float)
    n = code.size * spc
    fs = chip_rate * spc
    x_pm = np.asarray(x_pm if x_pm is not None else np.linspace(-650.0, 650.0, 64), float)
    elems, elems_r = [], []
    for s in sensors:
        L, k0, _, _ = calibrate(s["g"], s["R"])
        elems.append(dict(z=s["z"], lam_b=LAM0 + s["det"] * 1e-12, length=L, kappa0=k0, apod=GRATINGS[s["g"]].get("apod", "blackman")))
    for r in refs:
        L, k0, _, _ = calibrate(r["g"], r["R"])
        elems_r.append(dict(z=r["z"], lam_b=LAM0 + r["det"] * 1e-12, length=L, kappa0=k0, apod=GRATINGS[r["g"]].get("apod", "blackman")))
    I = drive_current(code, spc, bias, amp, rise_frac)
    fk = np.fft.fftfreq(n, 1.0 / fs)
    b_, a_ = bessel(4, 0.75 * chip_rate, fs=fs, norm="mag")
    det = np.zeros((len(seeds), x_pm.size, n), np.float32)
    for si, seed in enumerate(seeds):
        key = (seed, bias, amp, rise_frac, fs, code.size, int(code.sum()), str(laser_params), str(thermal))
        if key not in _LCACHE:
            _LCACHE[key] = laser_record(I, fs, seed=seed, params=laser_params, thermal=thermal)
        E_L, f_c = _LCACHE[key]
        rng = np.random.default_rng(seed * 7919)
        for m, x in enumerate(x_pm):
            lam_m = LAM0 + x * 1e-12
            nu = C0 / lam_m + f_c + fk
            lam = C0 / nu
            refl = array_reflection if ghosts else array_reflection_direct
            r = refl(lam, elems)
            if elems_r:
                r = r + refl(lam, elems_r)
            E_in = E_L * np.exp(1j * phase_walk(n, fs, linewidth_hz, rng))
            E_out = np.fft.ifft(np.fft.fft(E_in) * r)
            i_pd = photocurrent(E_out, fs, det_noise=det_noise, rng=rng)
            det[si, m] = lfilter(b_, a_, np.tile(i_pd, 2))[n:]      # causal filter, periodic warm-up
            if verbose and (m % 16 == 0 or m == x_pm.size - 1):
                print("   seed %d step %d/%d" % (seed, m + 1, x_pm.size), flush=True)
    meta = dict(z=np.array([s["z"] for s in sensors]), det_pm=np.array([s["det"] for s in sensors]), R=np.array([s["R"] for s in sensors]),
                gtype=np.array([s["g"] for s in sensors]), zr=np.array([r["z"] for r in refs]), detr_pm=np.array([r["det"] for r in refs]),
                Rr=np.array([r["R"] for r in refs]), code=code.astype(int), chip_rate=chip_rate, n_group=N_GROUP, lam0=LAM0,
                x_pm=x_pm, seeds=np.array(seeds), spc=spc, n_chips=code.size, linewidth_hz=linewidth_hz)
    if out:
        np.savez_compressed(out, det=det, **meta)
    return det, meta


def _mls(nbits, taps=None):
    taps = taps or {7: (7, 6), 9: (9, 5)}[nbits]
    n = 2 ** nbits - 1; reg = [1] * nbits; out = []
    for _ in range(n):
        out.append(reg[-1]); fb = 0
        for t in taps:
            fb ^= reg[t - 1]
        reg = [fb] + reg[:-1]
    return np.array(out)


if __name__ == "__main__":
    import time
    for g, R in (("bl250", 0.10), ("bl100", 0.01), ("uni40", 0.10)):
        L, k0, pk, fw = calibrate(g, R)
        print("%s R=%.2f: length %.2f mm, kappa0 %.0f 1/m, peak %.4f, FWHM %.0f pm" % (g, R, L * 1e3, k0, pk, fw))
    t0 = time.time()
    I = drive_current(_mls(7), 64, 0.033, 0.018)
    E, f_c = laser_record(I, 25e6 * 64, seed=1)
    P = np.abs(E) ** 2
    print("laser: %.1f s, on/off power ratio %.2f, mean |E|^2 %.3g" % (time.time() - t0, P[I > 0.033].mean() / P[I < 0.033].mean(), P.mean()))
    t0 = time.time()
    det, meta = run_array([dict(z=4.0, det=0.0, R=0.10, g="bl250")], [], x_pm=np.linspace(-650, 650, 64), seeds=(1,))
    print("single grating sweep: %.1f s, record max %.3e" % (time.time() - t0, det.max()))
