"""s11_sensitivity.py - how robust is the reference-corrected residual to the two
modelling choices it depends on: the assumed spread of sensor lineshape asymmetry,
and the assumed lineshape MODEL itself.

The headline residual (~3 pm after correction) is only meaningful if it does not
swing wildly when we change plausible assumptions. Two are tested:

  (a) sensor lineshape-asymmetry spread. Sensors are drawn with asymmetry in
      [-A, +A]; we sweep A over 0..0.4 and report the reference-corrected residual
      RMS. A larger spread means more slope mismatch between sensor and reference,
      so the residual should grow roughly linearly - the slope [pm per unit spread]
      is the sensitivity, and it should be modest.

  (b) lineshape model. The default is tanh^2 (strong-grating, coupled-mode). Weak
      gratings are closer to sinc^2 (Fourier transform of a uniform grating), which
      has visible sidelobes. We implement a simple normalized sinc^2 and compare,
      for the SAME chirp, (i) the raw FM-AM apparent error and (ii) the Gaussian-fit
      bias each model produces, plus the reference-corrected residual. The point is
      to check the conclusion is not an artifact of the tanh^2 choice.

sinc^2 is defined here locally (common.py is not modified). Its FWHM is matched to
FBG_FWHM_GHZ so the two models are compared at equal linewidth.
"""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import warnings; warnings.filterwarnings('ignore'); import common as C
np.random.seed(11)

PM = C.PM_PER_GHZ
F  = C.FBG_FWHM_GHZ                 # 250 pm

# ---------------------------------------------------------------------------
# weak-grating sinc^2 lineshape (visible sidelobes), normalized, FWHM-matched.
# sinc^2(x) = (sin(pi x)/(pi x))^2 has FWHM ~ 0.8859 in x; scale so the main lobe
# FWHM equals `fwhm`. Optional mild asymmetry (edge tilt) mirrors fbg_tanh's n_side.
# ---------------------------------------------------------------------------
_SINC2_FWHM_X = 0.8858929413  # FWHM of sinc^2 in units of x (numerically)
def fbg_sinc2(nu, nu_b, fwhm, n_side=0.0):
    x = (nu - nu_b) / (fwhm / _SINC2_FWHM_X)
    base = np.sinc(x)**2                       # np.sinc(x) = sin(pi x)/(pi x)
    if n_side != 0.0:
        base = base * (1.0 + n_side * np.tanh(x))
        base = np.clip(base, 0.0, None)
    return base / base.max()

def fmam_readout_shape(nu_grid, nu_b, delta, asym, shapefun, mean_off, skew=1.2):
    """fmam_readout equivalent but with an arbitrary lineshape function.
    Mirrors common.fmam_readout (same chirp_kernel, same mean_off handling)."""
    d, p = C.chirp_kernel(delta, skew=skew)
    d = d + mean_off
    out = np.zeros_like(nu_grid)
    for di, pi in zip(d, p):
        out += pi * shapefun(nu_grid + di, nu_b, F, n_side=asym)
    return out

def app_err_shape(nu_b, delta, asym, shapefun, off=True, skew=1.2):
    """Apparent Bragg error [pm] for a given lineshape model."""
    g = np.linspace(nu_b - 8*F, nu_b + 8*F, 1600)
    mo = 0.30*delta if off else 0.0
    spec = fmam_readout_shape(g, nu_b, delta, asym, shapefun, mean_off=mo, skew=skew)
    return (C.gauss_fit_peak(g, spec) - nu_b) * PM

def calibrate(se, sx, re, rx, order=1):
    order = min(order, len(rx)-1)
    return se - np.polyval(np.polyfit(rx, re, order), sx)

# ===========================================================================
# (a) reference-corrected residual vs sensor asymmetry SPREAD (0..0.4)
# ===========================================================================
delta = C.MEAS_CHIRP_GHZ
nsen = 12
ref_nu = np.array([-18.0, +18.0]); ref_as = np.array([+0.01, -0.01])
rng = np.random.default_rng(11)
sen_nu = np.sort(rng.uniform(-25.0, 25.0, nsen))

spreads = np.linspace(0.0, 0.4, 9)
res_tanh = []; res_sinc = []
for A in spreads:
    # draw sensor asymmetries within [-A,+A] (random sign pattern per spread level)
    signs = np.where(rng.uniform(-1, 1, nsen) >= 0, 1.0, -1.0)
    sen_as = A * signs
    # tanh^2 model
    se_t = np.array([app_err_shape(sen_nu[i], delta, sen_as[i], C.fbg_tanh) for i in range(nsen)])
    re_t = np.array([app_err_shape(ref_nu[j], delta, ref_as[j], C.fbg_tanh) for j in range(2)])
    res_tanh.append(np.sqrt(np.mean(calibrate(se_t, sen_nu, re_t, ref_nu)**2)))
    # sinc^2 model
    se_s = np.array([app_err_shape(sen_nu[i], delta, sen_as[i], fbg_sinc2) for i in range(nsen)])
    re_s = np.array([app_err_shape(ref_nu[j], delta, ref_as[j], fbg_sinc2) for j in range(2)])
    res_sinc.append(np.sqrt(np.mean(calibrate(se_s, sen_nu, re_s, ref_nu)**2)))
res_tanh = np.array(res_tanh); res_sinc = np.array(res_sinc)
slope_tanh = np.polyfit(spreads, res_tanh, 1)[0]
slope_sinc = np.polyfit(spreads, res_sinc, 1)[0]

# ===========================================================================
# (b) lineshape model comparison: tanh^2 vs sinc^2 at matched FWHM.
#   (i)  raw FM-AM apparent error vs chirp (symmetric grating, asym=0)
#   (ii) Gaussian-fit bias on the NOISELESS, UNCHIRPED lineshape (pure fit bias:
#        a Gaussian does not perfectly match either tail, most of all sinc^2's
#        sidelobes -> a static bias independent of chirp)
# ===========================================================================
deltas = np.linspace(0.3, 14.0, 20)   # chirp excursion [GHz]
raw_tanh = np.array([abs(app_err_shape(0.0, d, 0.0, C.fbg_tanh, off=False)) for d in deltas])
raw_sinc = np.array([abs(app_err_shape(0.0, d, 0.0, fbg_sinc2, off=False)) for d in deltas])

def gauss_fit_bias(shapefun):
    """Gaussian-fit peak offset [pm] on the noiseless, unchirped lineshape."""
    g = np.linspace(-8*F, 8*F, 3201)
    y = shapefun(g, 0.0, F, n_side=0.0)
    return C.gauss_fit_peak(g, y) * PM      # should be ~0 by symmetry; residual = fit bias
bias_tanh = gauss_fit_bias(C.fbg_tanh)
bias_sinc = gauss_fit_bias(fbg_sinc2)

# widths a Gaussian fit assigns to each (diagnostic: how badly Gauss matches sinc^2)
def gauss_fit_sigma(shapefun):
    g = np.linspace(-8*F, 8*F, 3201); y = shapefun(g, 0.0, F, n_side=0.0)
    i0 = int(np.argmax(y)); W = int(len(g)*0.14); lo=max(0,i0-W); hi=min(len(g),i0+W+1)
    xs=g[lo:hi]; ys=y[lo:hi]
    try:
        p,_=curve_fit(lambda x,a,mu,s,c: a*np.exp(-0.5*((x-mu)/s)**2)+c, xs, ys,
                      p0=[ys.max(), 0.0, F/2.0, 0.0], maxfev=20000)
        return abs(p[2])*PM
    except Exception:
        return np.nan
sig_tanh = gauss_fit_sigma(C.fbg_tanh); sig_sinc = gauss_fit_sigma(fbg_sinc2)

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(13.4, 3.7))
# (a) residual vs asymmetry spread, both models
ax[0].plot(spreads, res_tanh, 'o-', color='#2980b9', label=f'tanh$^2$ ({slope_tanh:.1f} pm/unit)')
ax[0].plot(spreads, res_sinc, 's-', color='#c0392b', label=f'sinc$^2$ ({slope_sinc:.1f} pm/unit)')
ax[0].set_xlabel('sensor asymmetry spread A'); ax[0].set_ylabel('corrected residual RMS [pm]')
ax[0].set_title('(a) Residual vs asymmetry spread'); ax[0].legend(fontsize=7)
# (b) the two lineshapes (show sidelobes)
gg = np.linspace(-4*F, 4*F, 1600)
ax[1].plot(gg*PM/1000, C.fbg_tanh(gg, 0.0, F), color='#2980b9', label='tanh$^2$')
ax[1].plot(gg*PM/1000, fbg_sinc2(gg, 0.0, F), color='#c0392b', label='sinc$^2$ (weak grating)')
ax[1].set_xlabel('wavelength offset [nm]'); ax[1].set_ylabel('norm. reflectivity')
ax[1].set_title('(b) Lineshape models (FWHM-matched)'); ax[1].legend(fontsize=7)
ax[1].set_ylim(-0.02, 1.05)
# (c) raw FM-AM apparent error vs chirp, both models
ax[2].plot(deltas/F, raw_tanh, 'o-', color='#2980b9', label='tanh$^2$')
ax[2].plot(deltas/F, raw_sinc, 's-', color='#c0392b', label='sinc$^2$')
ax[2].axvspan(C.CHIRP_REALISTIC_GHZ[0]/F, C.CHIRP_REALISTIC_GHZ[1]/F, color='0.9')
ax[2].set_xlabel('chirp excursion / FWHM'); ax[2].set_ylabel('|raw FM-AM error| [pm]')
ax[2].set_title('(c) Raw FM-AM error vs chirp per model'); ax[2].legend(fontsize=7)
plt.tight_layout(); plt.savefig('figs/fig_s11_sensitivity.png', dpi=140)

# ---------------------------------------------------------------------------
# printed table
# ---------------------------------------------------------------------------
print(f"chirp (for panels a,b) = {delta:.1f} GHz; FWHM = {C.FBG_FWHM_PM:.0f} pm; 2 references")
print("--- (a) reference-corrected residual RMS vs sensor asymmetry spread ---")
print("   spread A   tanh^2 [pm]   sinc^2 [pm]")
for A, rt, rs in zip(spreads, res_tanh, res_sinc):
    print(f"   {A:6.2f}     {rt:8.2f}     {rs:8.2f}")
print(f"  slope: tanh^2 = {slope_tanh:.2f} pm/unit spread, sinc^2 = {slope_sinc:.2f} pm/unit spread")
print("--- (b) lineshape-model comparison at matched FWHM ---")
print(f"  Gaussian-fit bias on the symmetric, unchirped lineshape:")
print(f"     tanh^2 = {bias_tanh:+.3f} pm   |   sinc^2 = {bias_sinc:+.3f} pm")
print(f"  Gaussian-fit sigma assigned (main-lobe width proxy):")
print(f"     tanh^2 = {sig_tanh:.1f} pm     |   sinc^2 = {sig_sinc:.1f} pm")
i_meas = int(np.argmin(abs(deltas - C.MEAS_CHIRP_GHZ)))
print(f"  raw FM-AM error @measured chirp {C.MEAS_CHIRP_GHZ:.1f} GHz:")
print(f"     tanh^2 = {raw_tanh[i_meas]:.1f} pm   |   sinc^2 = {raw_sinc[i_meas]:.1f} pm")
print("saved figs/fig_s11_sensitivity.png")
