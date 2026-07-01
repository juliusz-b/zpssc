"""s8_noise.py - detector-noise floor of the Bragg-position readout, and where
it sits relative to the systematic FM-AM residual.

Every result so far (s3-s7) is a SYSTEMATIC error: the chirp/lineshape residual
that survives reference correction (~3 pm). That residual is deterministic - it
does not average down with more trials. Detector noise is the opposite: a random
per-shot jitter that DOES average down. This script asks a practical question:
at what electrical SNR does the random centroid jitter fall below the systematic
residual, so that noise stops being the limiting term?

Model. Take a single FBGS grating (tanh^2, FWHM = common.FBG_FWHM_GHZ = 250 pm,
R = 0.10). The despread spectral readout is the reflected power vs sweep position.
We add additive Gaussian noise referred to that readout at a controlled ELECTRICAL
SNR (SNR = signal_power / noise_power, in dB, defined on the peak amplitude), fit
the Bragg position with common.gauss_fit_peak over many trials, and report the RMS
centroid jitter [pm] vs SNR. Lower R means a weaker return, so at fixed SNR the
jitter grows as the peak amplitude shrinks - shown in the second panel.

Note on framing: the ~3 pm line is the reference-corrected FM-AM residual from
s3/s6 (the genuine edge distortion, not the removable common offset). It is drawn
as a horizontal marker so the reader sees the SNR at which noise dominates.
"""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore'); import common as C
np.random.seed(8)

PM = C.PM_PER_GHZ
F  = C.FBG_FWHM_GHZ                 # ~31.25 GHz (= 250 pm)
R  = C.FBG_R                        # 0.10
NU_B = 0.0                          # single grating at band center
SYS_RESIDUAL_PM = 3.0              # reference-corrected FM-AM residual (s3/s6 headline)

# Sweep grid: several FWHM either side, sampled finely enough that the fit is not
# grid-limited (grid step ~ F/75 -> sub-pm quantization of the fitted peak). The
# noise, not the grid, dominates the jitter, so a moderate grid is sufficient.
nu = np.linspace(NU_B - 6*F, NU_B + 6*F, 900)

def readout(nu_b=NU_B, refl=R, asym=0.0):
    """Noiseless despread spectral readout (reflected power), peak amplitude = refl."""
    return refl * C.fbg_tanh(nu, nu_b, F, n_side=asym)

def jitter_sigma(sigma, refl=R, ntrials=400, asym=0.0, rng=None):
    """RMS Bragg-position jitter [pm] for an ABSOLUTE additive-noise std `sigma`
    referred to the readout (electronics/thermal floor, independent of signal)."""
    if rng is None: rng = np.random.default_rng(0)
    clean = readout(refl=refl, asym=asym)
    est = np.empty(ntrials)
    for k in range(ntrials):
        noisy = clean + rng.normal(0.0, sigma, size=clean.shape)
        est[k] = C.gauss_fit_peak(nu, noisy)
    ref_pos = C.gauss_fit_peak(nu, clean)     # remove any static grid/fit bias
    return np.std((est - ref_pos) * PM)

def snr_jitter(snr_db, refl=R, ntrials=400, asym=0.0, rng=None):
    """Jitter at an electrical SNR defined on the readout PEAK amplitude of this
    grating: sigma = peak / 10**(SNR/20). Used for the SNR sweep (panel a)."""
    peak = readout(refl=refl, asym=asym).max()
    sigma = peak / (10.0**(snr_db/20.0))
    return jitter_sigma(sigma, refl=refl, ntrials=ntrials, asym=asym, rng=rng)

# ---------------------------------------------------------------------------
# (a) jitter vs SNR at the procured reflectivity R = 0.10
# ---------------------------------------------------------------------------
snr_db = np.arange(10.0, 41.0, 2.5)
rng = np.random.default_rng(8)
NTRIAL = 300
jit_snr = np.array([snr_jitter(s, refl=R, ntrials=NTRIAL, rng=rng) for s in snr_db])

# SNR at which random jitter crosses below the systematic residual
snr_cross = float(np.interp(SYS_RESIDUAL_PM, jit_snr[::-1], snr_db[::-1]))  # monotone decr.

# ---------------------------------------------------------------------------
# (b) jitter vs reflectivity at a fixed SNR (weaker return -> more jitter)
# ---------------------------------------------------------------------------
SNR_FIX = 25.0
Rs = np.linspace(0.02, 0.30, 11)
rng2 = np.random.default_rng(80)
# Absolute noise floor: fix sigma from SNR_FIX at the procured R=0.10, then keep it
# CONSTANT while R varies. A weaker grating (lower R) has a lower peak, so its
# effective SNR degrades and the jitter grows ~ 1/R.
peak_ref = readout(refl=R).max()
sigma_floor = peak_ref / (10.0**(SNR_FIX/20.0))
jit_R = np.array([jitter_sigma(sigma_floor, refl=rr, ntrials=NTRIAL, rng=rng2) for rr in Rs])

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(10.6, 3.8))
ax[0].semilogy(snr_db, jit_snr, 'o-', label='RMS centroid jitter (Gauss fit)')
ax[0].axhline(SYS_RESIDUAL_PM, color='#c0392b', ls='--', lw=1.2,
              label='systematic FM-AM residual (~%.0f pm)' % SYS_RESIDUAL_PM)
ax[0].axvline(snr_cross, color='0.5', ls=':', lw=1.0)
ax[0].text(snr_cross, jit_snr.max()*0.55, ' noise = residual\n @%.0f dB' % snr_cross,
           fontsize=6.5, ha='left')
ax[0].set_xlabel('electrical SNR [dB]'); ax[0].set_ylabel('RMS Bragg jitter [pm]')
ax[0].set_title('(a) Position jitter vs SNR (R = 0.10, FWHM 250 pm)')
ax[0].legend(fontsize=7); ax[0].grid(True, which='both', alpha=0.25)
ax[1].plot(Rs*100, jit_R, 's-', color='#2980b9', label='jitter @SNR = %.0f dB' % SNR_FIX)
ax[1].axhline(SYS_RESIDUAL_PM, color='#c0392b', ls='--', lw=1.2,
              label='systematic residual (~%.0f pm)' % SYS_RESIDUAL_PM)
ax[1].axvline(R*100, color='0.5', ls=':', lw=1.0)
ax[1].text(R*100, jit_R.max()*0.9, ' procured\n R = %.2f' % R, fontsize=6.5, ha='left')
ax[1].set_xlabel('grating reflectivity R [%]'); ax[1].set_ylabel('RMS Bragg jitter [pm]')
ax[1].set_title('(b) Jitter vs reflectivity (SNR = %.0f dB)' % SNR_FIX)
ax[1].legend(fontsize=7); ax[1].grid(True, alpha=0.25)
plt.tight_layout(); plt.savefig('figs/fig_s8_noise.png', dpi=140)

# ---------------------------------------------------------------------------
# printed tables
# ---------------------------------------------------------------------------
print("single FBG: nu_b=%.0f GHz, FWHM=%.0f pm (=%.1f GHz), R=%s" % (NU_B, C.FBG_FWHM_PM, F, R))
print("systematic FM-AM residual reference line: %.1f pm" % SYS_RESIDUAL_PM)
print("--- (a) RMS Bragg jitter vs electrical SNR (R = 0.10) ---")
print("   SNR[dB]   jitter[pm]   vs residual")
for s, j in zip(snr_db, jit_snr):
    tag = "noise-limited" if j > SYS_RESIDUAL_PM else "residual-limited"
    print("   %6.0f   %9.2f   %s" % (s, j, tag))
print("  => random jitter drops below the ~%.0f pm systematic residual at SNR ~ %.0f dB;"
      % (SYS_RESIDUAL_PM, snr_cross))
print("     above that, the systematic (not noise) sets the floor.")
print("--- (b) RMS Bragg jitter vs reflectivity (SNR = %.0f dB) ---" % SNR_FIX)
print("     R[%]   jitter[pm]")
for rr, j in zip(Rs, jit_R):
    print("   %6.1f   %9.2f" % (rr*100, j))
jit_at_R = float(np.interp(R, Rs, jit_R))
print("  => at the procured R=%.2f, jitter ~ %.2f pm at %.0f dB; halving R roughly doubles the jitter."
      % (R, jit_at_R, SNR_FIX))
print("saved figs/fig_s8_noise.png")
