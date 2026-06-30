"""
common.py - shared physics for the simulation study: self-calibrating CDM
interrogation of a fiber Bragg grating array with a directly-modulated, swept
VCSEL. Project PN/01/0321/2022 (ZPSSC).

Units: optical frequency in GHz relative to band center; wavelength shifts in pm,
conversion 1 GHz ~ 8 pm at 1550 nm.
"""
import numpy as np
from scipy.signal import max_len_seq
from scipy.optimize import curve_fit
import warnings; warnings.filterwarnings("ignore")

GHZ_PER_PM = 0.125            # 1 pm -> 0.125 GHz @1550 nm  (dnu = c*dlam/lam^2)
PM_PER_GHZ = 1.0 / GHZ_PER_PM # ~8.0 pm per GHz

# ----------------------------------------------------------------------------
# Procured hardware (sets experiment-realistic defaults so the simulation maps
# onto the actual bench): FBGS DTG-A3A4 gratings, order S-2026-0066, and the
# 3-channel Peltier stage. 3 gratings used on the bench, 15 procured. All
# gratings share the SAME nominal wavelength -> they overlap spectrally and are
# separated only by code/delay (the core CDM claim). A small wavelength band is
# created by tuning each grating's temperature (Peltier).
# ----------------------------------------------------------------------------
LAMBDA0_NM         = 1545.0    # nominal Bragg wavelength of every grating (identical)
FBG_FWHM_PM        = 250.0     # FBGS DTG declared linewidth (FWHM); confirm on test sheet
FBG_FWHM_GHZ       = FBG_FWHM_PM * GHZ_PER_PM   # ~31.25 GHz
FBG_R              = 0.10      # reflectivity per grating (config A0LF, declared 10%)
N_GRATINGS_BENCH   = 3         # gratings on the bench (one Peltier stage each)
N_GRATINGS_AVAIL   = 15        # gratings procured (order S-2026-0066)
TEMP_COEF_PM_PER_C = 10.0      # Bragg temperature sensitivity (~10 pm/C)
TEMP_RANGE_C       = (15.0, 60.0)   # Peltier set-point range (to confirm)
TUNING_RANGE_PM    = (TEMP_RANGE_C[1] - TEMP_RANGE_C[0]) * TEMP_COEF_PM_PER_C  # ~450 pm
TEMP_STAB_C        = 0.1       # PID stability -> ~1 pm wavelength stability

def temp_to_dnu_ghz(dT_celsius):
    """Bragg detuning (GHz) for a temperature change dT (Celsius) at ~10 pm/C."""
    return dT_celsius * TEMP_COEF_PM_PER_C * GHZ_PER_PM

# ----------------------------------------------------------------------------
# Spreading sequences (m-sequences, Gold, Kasami) - decimation construction
# ----------------------------------------------------------------------------
def _mls01(nbits):
    seq, _ = max_len_seq(nbits)      # 0/1, length N=2^n-1
    return seq.astype(int)

def _decimate01(seq, q):
    N = len(seq)
    idx = (q * np.arange(N)) % N
    return seq[idx]

def _to_pm1(seq01):
    return 1.0 - 2.0 * seq01         # 0->+1, 1->-1  (bipolar)

def gold_set(nbits):
    """Gold family for odd n (preferred pair via decimation q=2^((n+1)/2)+1)."""
    N = 2**nbits - 1
    u = _mls01(nbits)
    q = 2**((nbits + 1)//2) + 1
    v = _decimate01(u, q)
    codes = [u.copy(), v.copy()]
    for k in range(N):
        codes.append((u ^ np.roll(v, k)))
    return np.array([_to_pm1(c) for c in codes])

def kasami_small(nbits):
    """Small Kasami set for even n (decimation q=2^(n/2)+1)."""
    assert nbits % 2 == 0
    N = 2**nbits - 1
    u = _mls01(nbits)
    q = 2**(nbits//2) + 1
    w = _decimate01(u, q)            # period 2^(n/2)-1, tiled to N
    codes = [u.copy()]
    for k in range(2**(nbits//2) - 1):
        codes.append(u ^ np.roll(w, k))
    return np.array([_to_pm1(c) for c in codes])

def periodic_xcorr(a, b):
    """Cyclic cross-correlation (normalized by length)."""
    n = len(a)
    A = np.fft.fft(a); B = np.fft.fft(b)
    return np.fft.ifft(A * np.conj(B)).real / n

def code_corr_stats(codes):
    """Peak auto-correlation side lobe and peak cross-correlation (norm. to N)."""
    M, N = codes.shape
    auto_side = 0.0; cross = 0.0
    for i in range(M):
        c = periodic_xcorr(codes[i], codes[i])
        auto_side = max(auto_side, np.max(np.abs(c[1:])))   # excluding the zero-lag peak
        for j in range(i+1, M):
            cc = periodic_xcorr(codes[i], codes[j])
            cross = max(cross, np.max(np.abs(cc)))
    return auto_side, cross   # zero-lag peaks equal 1 (=N/N)

# ----------------------------------------------------------------------------
# Fiber Bragg grating spectra
# ----------------------------------------------------------------------------
def fbg_gauss(nu, nu_b, fwhm):
    sig = fwhm / 2.35482
    return np.exp(-0.5 * ((nu - nu_b) / sig)**2)

def fbg_tanh(nu, nu_b, fwhm, n_side=0.0):
    """tanh^2 (coupled-mode) approximation with optional asymmetry n_side.
    The asymmetry models the real, slightly non-symmetric grating lineshape."""
    sig = fwhm / 2.35482
    x = (nu - nu_b) / sig
    base = 1.0 / np.cosh(x)**2
    if n_side != 0.0:
        base = base * (1.0 + n_side * np.tanh(x))   # mild edge asymmetry
        base = np.clip(base, 0.0, None)
    return base / base.max()

# ----------------------------------------------------------------------------
# Chirp / FM-to-AM conversion
# ----------------------------------------------------------------------------
def chirp_kernel(delta, skew=0.0, npts=401, span=4.0):
    """Distribution of instantaneous frequency offsets during the 'on' code chips.
    delta - magnitude (std) of the chirp excursion [GHz]; skew - skewness (direct
    modulation gives an asymmetric transient-vs-adiabatic chirp)."""
    if delta <= 0:
        d = np.array([0.0]); p = np.array([1.0]); return d, p
    d = np.linspace(-span*delta, span*delta, npts)
    p = np.exp(-0.5 * (d/delta)**2)
    if skew != 0.0:
        from scipy.special import erf
        p = p * (1.0 + erf(skew * d / (np.sqrt(2)*delta)))   # skewed distribution
    p = np.clip(p, 0, None); p /= p.sum()
    return d, p

def fmam_readout(nu_grid, nu_b, fwhm, delta, mean_off=0.0, skew=0.0,
                 shape='tanh', asym=0.0):
    """Effective despread spectral readout of a grating with FM->AM blur.
    The chirp makes the receiver, at sweep position nu_s, average the grating
    spectrum over the distribution of instantaneous frequencies (chirp kernel)."""
    d, p = chirp_kernel(delta, skew=skew)
    d = d + mean_off
    if shape == 'gauss':
        spec = lambda x: fbg_gauss(x, nu_b, fwhm)
    else:
        spec = lambda x: fbg_tanh(x, nu_b, fwhm, n_side=asym)
    out = np.zeros_like(nu_grid)
    for di, pi in zip(d, p):
        out += pi * spec(nu_grid + di)
    return out

# ----------------------------------------------------------------------------
# Peak-position estimators
# ----------------------------------------------------------------------------
def _gauss(x, a, mu, sig, c):
    return a * np.exp(-0.5*((x-mu)/sig)**2) + c

def centroid(x, y, thr=0.3):
    ym = y.max(); m = y >= thr*ym
    w = y[m] - thr*ym
    return np.sum(x[m]*w)/np.sum(w)

def gauss_fit_peak(x, y, thr=0.3, win_frac=0.14):
    """Gaussian fit in a local window around the peak (robust to side lobes/cross-talk)."""
    x=np.asarray(x,float); y=np.asarray(y,float)
    i0=int(np.argmax(y)); W=max(4,int(len(x)*win_frac))
    lo=max(0,i0-W); hi=min(len(x),i0+W+1)
    xs=x[lo:hi]; ys=y[lo:hi]
    ym=ys.max(); base=float(np.min(ys))
    m=ys>=base+thr*(ym-base)
    if m.sum()<4: return centroid(xs,ys,thr)
    p0=[ym-base, xs[int(np.argmax(ys))], (xs[-1]-xs[0])/6.0, base]
    try:
        lob=[0.0,xs[0],1e-6,-np.inf]; hib=[np.inf,xs[-1],xs[-1]-xs[0],np.inf]
        popt,_=curve_fit(_gauss,xs[m],ys[m],p0=p0,bounds=(lob,hib),maxfev=20000)
        return float(popt[1])
    except Exception:
        return centroid(xs,ys,thr)

# ----------------------------------------------------------------------------
# Nonlinear VCSEL sweep and Mach-Zehnder interferometer (k-clock)
# ----------------------------------------------------------------------------
def vcsel_sweep(v, span_ghz, quad=0.30, cube=0.0):
    """nu(V) ~ quadratic-cubic function of the HCG-VCSEL tuning voltage
    (BW10: delta lambda ~ quadratic in tuning voltage). v in [0,1]."""
    base = v + quad*(v**2 - v) + cube*(v**3 - v)
    base = (base - base.min())/(base.max()-base.min())
    return base * span_ghz

def mzi_kclock_estimate(nu_true, fsr_ghz):
    """MZI frequency ruler: returns 'ticks' at equal increments of nu (fringe
    zeros). An ideal ruler of the source frequency (including the common chirp),
    but blind to the grating-edge FM->AM (a detection-side effect)."""
    n0 = nu_true.min(); n1 = nu_true.max()
    ticks = np.arange(n0, n1, fsr_ghz)
    return ticks
