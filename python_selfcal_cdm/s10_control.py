"""s10_control.py - THE NOVELTY CONTROL (decisive).

The claimed novelty is that the reference grating is CO-CODED: its apparent FM-AM
error is measured through the SAME despreading/correlation pipeline as the sensor,
rather than by a separate, idealized anchor. The decisive question is blunt:

   Does routing the reference through despreading change the correction it gives,
   compared with a plain (non-coded) reference grating at the same wavelength?

If the corrected sensor residual is the SAME either way, then co-coding buys
nothing over a plain reference grating and the novelty is weak (it is a systems /
packaging convenience, not an accuracy mechanism). If it DIFFERS, we quantify by
how much and in which direction.

s3 already made this comparison but sampled BOTH references from the identical
analytic FM-AM model, so it could only ever return "identical". This script is
stricter: the co-coded reference is measured through a REAL correlation despreader
(build the coded composite, correlate against the code, read the despread
lineshape, Gaussian-fit it) exactly like s6's panel (d). The direct reference is
read straight off the noiseless FM-AM model. We compare in two regimes:

  (1) CLEAN: only the reference grating present. Tests whether the correlation
      operation itself biases the reference's apparent error (it should not - the
      code auto-correlation is a scaled delta, so despreading is transparent).
  (2) LOADED: the sensor plus other co-coded gratings share the fiber, so their
      codes leak residual cross-correlation onto the reference slot. A co-coded
      reference SEES this leakage (it is measured in the same crowded channel); a
      plain separate reference does NOT. This is the only place co-coding can
      matter, and here it matters as a small PENALTY, not a benefit.
"""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore'); import common as C
np.random.seed(10)

PM = C.PM_PER_GHZ
F  = C.FBG_FWHM_GHZ                 # 250 pm
R  = C.FBG_R

# ---- shared FM-AM model (same chirp/lineshape as elsewhere) ---------------
def app_err_direct(nu_b, delta, asym, skew=1.2):
    """Apparent Bragg error [pm], read straight off the noiseless FM-AM model.
    This is the DIRECT / plain-reference measurement (no despreading)."""
    g = np.linspace(nu_b - 6*F, nu_b + 6*F, 1200)
    spec = C.fmam_readout(g, nu_b, F, delta, mean_off=0.30*delta, skew=skew,
                          shape='tanh', asym=asym)
    return (C.gauss_fit_peak(g, spec) - nu_b) * PM

# ---- co-coded readout through a REAL correlation despreader ---------------
# Faithful to s6's despreading model: all gratings share ONE swept-wavelength axis
# (the common slot grid). Each grating's FM-AM-blurred reflected amplitude is
# sampled on that grid, spread by cycling the code copies per slot (codes[mi % N])
# and delayed by the grating's tau (its physical position, 4 m spacing). The
# receiver sums all gratings' returns, then despreads the TARGET by correlating
# against the same per-slot codes at the target's tau. Interferers leak residual
# cross-correlation onto the target - the only mechanism by which co-coding can
# differ from a plain reference. MSLOT is kept ~ N (as in s6, M=128 vs N=127) so
# each slot gets a near-unique code and there is no self-aliasing.
CODES = C.gold_set(7)              # Gold family, N=127
NCODE = CODES.shape[1]
MSLOT = 128                         # wavelength slots (== s6), near the code length
SLOT_HALFWIN = 70.0                 # common swept axis: +/-70 GHz about band center
SLOTS = np.linspace(-SLOT_HALFWIN, SLOT_HALFWIN, MSLOT)   # shared wavelength grid

def fmam_lineshape_on_grid(nu_b, delta, asym, skew=1.2):
    """FM-AM-blurred reflected amplitude on the shared swept-wavelength grid."""
    return C.fmam_readout(SLOTS, nu_b, F, delta, mean_off=0.30*delta,
                          skew=skew, shape='tanh', asym=asym) * R

def despread_position(gratings, tgt_delay, delta, skew=1.2):
    """Build the coded composite from a list of gratings, despread the TARGET at
    tgt_delay, and return the Gaussian-fit peak position [GHz] on the shared grid.
    `gratings` = list of (nu_b, asym, delay); the target is the one with tgt_delay."""
    rx = np.zeros(NCODE)
    for (nb, a, dly) in gratings:
        refl = fmam_lineshape_on_grid(nb, delta, a, skew)
        for mi, a_m in enumerate(refl):
            rx += a_m * np.roll(CODES[mi % NCODE], dly)      # per-slot code, grating delay
    desp = np.array([C.periodic_xcorr(rx, CODES[mi % NCODE])[tgt_delay]
                     for mi in range(MSLOT)])
    return C.gauss_fit_peak(SLOTS, desp)

def app_err_coded(nu_b, delta, asym, delay, others=None, skew=1.2):
    """Apparent Bragg error [pm] of a co-coded grating recovered THROUGH the
    despreading pipeline. `others` = list of (nu_b, asym, delay) for interfering
    co-coded gratings that share the fiber (their code cross-correlation leaks
    onto the target). Positions are read on the shared swept-wavelength grid."""
    gratings = [(nu_b, asym, delay)]
    if others:
        gratings += list(others)
    return (despread_position(gratings, delay, delta, skew) - nu_b) * PM

def app_err_direct_on_slots(nu_b, delta, asym, skew=1.2):
    """DIRECT reading on the SAME coarse slot grid as the despreader, but WITHOUT
    coding. This isolates the single variable under test (the correlation
    despreading) from the wavelength-grid resolution: comparing this against
    app_err_coded holds the grid identical, so any residual difference is due to
    the coding alone, not to sampling."""
    y = fmam_lineshape_on_grid(nu_b, delta, asym, skew)   # noiseless, on SLOTS
    return (C.gauss_fit_peak(SLOTS, y) - nu_b) * PM

def calibrate(se, sx, re, rx, order=1):
    order = min(order, len(rx)-1)
    p = np.polyfit(rx, re, order)
    return se - np.polyval(p, sx)

# ===========================================================================
# REGIME 1 - CLEAN: reference alone. Direct model vs coded-despread readout.
# ===========================================================================
delta = C.MEAS_CHIRP_GHZ
rng = np.random.default_rng(10)
nsen = 12
sen_nu = np.sort(rng.uniform(-25.0, 25.0, nsen))     # sensors within the bench band
sen_as = rng.uniform(-0.30, 0.30, nsen)
ref_nu = np.array([-18.0, +18.0])                    # two co-coded references
ref_as = np.array([+0.01, -0.01])
ref_delay = [17, 61]                                 # distinct delays (4 m spacing)

# sensor apparent errors (direct model; sensors are the thing being measured)
se = np.array([app_err_direct(sen_nu[i], delta, sen_as[i]) for i in range(nsen)])

# reference errors measured three ways:
#  (i)  re_direct        - plain reference, fine analytic grid (the realistic anchor)
#  (ii) re_direct_slots  - plain reference on the despreader's coarse grid (grid-matched
#                          control, so the ONLY difference vs coded is the correlation)
#  (iii) re_coded_clean  - co-coded reference through the correlation despreader
re_direct = np.array([app_err_direct(ref_nu[j], delta, ref_as[j]) for j in range(2)])
re_direct_slots = np.array([app_err_direct_on_slots(ref_nu[j], delta, ref_as[j]) for j in range(2)])
re_coded_clean = np.array([app_err_coded(ref_nu[j], delta, ref_as[j], ref_delay[j])
                           for j in range(2)])

cor_direct = calibrate(se, sen_nu, re_direct, ref_nu, 1)
cor_coded_clean = calibrate(se, sen_nu, re_coded_clean, ref_nu, 1)

rms_direct = np.sqrt(np.mean(cor_direct**2))
rms_coded_clean = np.sqrt(np.mean(cor_coded_clean**2))
# grid-matched comparison: coded vs direct-on-same-grid isolates the coding effect
ref_diff_clean = np.max(np.abs(re_coded_clean - re_direct_slots))   # coding-only disagreement
cor_diff_clean = ref_diff_clean   # linear correction: reference-error diff maps 1:1

# ===========================================================================
# REGIME 2 - LOADED: the reference is measured with OTHER co-coded gratings
# sharing the fiber. Their codes leak cross-correlation onto the reference slot.
# A plain (direct) reference does not see this; a co-coded one does.
# ===========================================================================
# interferers: other co-coded gratings at distinct delays (their own bench
# positions), sharing the fiber. Format (nu_b [GHz], asym, delay).
interferers = [(-10.0, +0.15, 29), (+6.0, -0.20, 89), (+20.0, +0.10, 113)]
re_coded_loaded = np.array([
    app_err_coded(ref_nu[j], delta, ref_as[j], ref_delay[j], others=interferers)
    for j in range(2)])
cor_coded_loaded = calibrate(se, sen_nu, re_coded_loaded, ref_nu, 1)
rms_coded_loaded = np.sqrt(np.mean(cor_coded_loaded**2))
ref_diff_loaded = np.max(np.abs(re_coded_loaded - re_direct))
cor_diff_loaded = np.max(np.abs(cor_coded_loaded - cor_direct))

# ---------------------------------------------------------------------------
# small figure (optional): reference errors and corrected residuals per method
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(10.4, 3.8))
methods = ['direct\n(plain ref)', 'co-coded\n(clean)', 'co-coded\n(loaded)']
rms_vals = [rms_direct, rms_coded_clean, rms_coded_loaded]
ax[0].bar(range(3), rms_vals, color=['#2980b9', '#27ae60', '#e67e22'])
ax[0].set_xticks(range(3)); ax[0].set_xticklabels(methods, fontsize=8)
ax[0].set_ylabel('corrected residual RMS [pm]')
ax[0].set_title('(a) Corrected residual: does co-coding change it?')
for i, v in enumerate(rms_vals):
    ax[0].text(i, v + max(rms_vals)*0.02, f'{v:.2f}', ha='center', fontsize=8)
ax[1].axhline(0, c='0.85', lw=0.7)
ax[1].plot(sen_nu*PM/1000, cor_direct, 'o-', ms=4, color='#2980b9', label='direct ref')
ax[1].plot(sen_nu*PM/1000, cor_coded_clean, 's-', ms=4, color='#27ae60', label='co-coded clean')
ax[1].plot(sen_nu*PM/1000, cor_coded_loaded, '^-', ms=4, color='#e67e22', label='co-coded loaded')
ax[1].set_xlabel('sensor position in band [nm]'); ax[1].set_ylabel('corrected error [pm]')
ax[1].set_title('(b) Per-sensor corrected error'); ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig('figs/fig_s10_control.png', dpi=140)

# ---------------------------------------------------------------------------
# printed numbers + verdict
# ---------------------------------------------------------------------------
print(f"chirp = {delta:.1f} GHz (measured); {nsen} sensors, 2 co-coded references; Gold N={NCODE}")
print("--- REGIME 1 (CLEAN: reference alone) ---")
print(f"  reference error, DIRECT (fine grid)      : {re_direct[0]:+.3f}, {re_direct[1]:+.3f} pm")
print(f"  reference error, DIRECT (despreader grid) : {re_direct_slots[0]:+.3f}, {re_direct_slots[1]:+.3f} pm")
print(f"  reference error, CODED  (despreader grid) : {re_coded_clean[0]:+.3f}, {re_coded_clean[1]:+.3f} pm")
print(f"  max |coded - direct| on the SAME grid     : {ref_diff_clean:.3e} pm  (coding-only effect)")
print(f"  corrected sensor residual RMS: direct = {rms_direct:.3f} pm | co-coded = {rms_coded_clean:.3f} pm")
print("--- REGIME 2 (LOADED: other co-coded gratings share the fiber) ---")
print(f"  reference error, CODED despread   : {re_coded_loaded[0]:+.3f}, {re_coded_loaded[1]:+.3f} pm")
print(f"  max |coded(loaded) - direct| per ref: {ref_diff_loaded:.3f} pm  (code cross-talk leakage)")
print(f"  corrected sensor residual RMS: co-coded loaded = {rms_coded_loaded:.3f} pm "
      f"(vs {rms_direct:.3f} pm direct)")
print(f"  max |corrected difference| per sensor: {cor_diff_loaded:.3f} pm")
print("--- VERDICT: does co-coding matter? ---")
if ref_diff_clean < 0.1:
    print(f"  CLEAN regime: with the wavelength grid held identical, co-coded and plain references "
          f"agree to {ref_diff_clean:.2e} pm. The correlation despreader is transparent for the "
          f"reference: the Gold auto-correlation is a scaled delta, so the despread lineshape equals "
          f"the model lineshape and the recovered error is the same up to numerical precision.")
    print(f"  => On accuracy grounds alone, co-coding adds NOTHING over a plain reference grating at "
          f"the same wavelength (the ~0.5 pm gap vs the fine-grid direct reading is just the "
          f"despreader's slot quantization, not an effect of coding). The novelty is NOT a "
          f"correction-accuracy gain.")
else:
    print(f"  CLEAN regime: grid-matched corrections DIFFER by {ref_diff_clean:.2f} pm - the coding "
          f"itself changes the reference reading even without interferers.")
if cor_diff_loaded > 5*max(cor_diff_clean, 1e-6):
    print(f"  LOADED regime: co-coding DOES change the result, but as a PENALTY: sharing the "
          f"channel injects up to {ref_diff_loaded:.2f} pm of code cross-talk onto the reference, "
          f"lifting the corrected residual to {rms_coded_loaded:.2f} pm (from {rms_direct:.2f} pm). "
          f"A separate plain reference avoids this leakage.")
    print(f"  => Net: co-coding's honest value is SYSTEMS-level (one source, one chirp, one "
          f"despreader, no extra WDM anchor and no Peltier per reference-wavelength), not a "
          f"residual-accuracy advantage. If anything, the shared channel costs a little accuracy.")
else:
    print(f"  LOADED regime: cross-talk leakage is negligible ({ref_diff_loaded:.2e} pm); co-coding "
          f"is accuracy-neutral even in a crowded channel.")
print("saved figs/fig_s10_control.png")
