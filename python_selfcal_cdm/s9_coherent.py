"""s9_coherent.py - source coherence length, and the extra error that coherent
beating could add if two co-located gratings' reflections combine in FIELD rather
than in POWER.

Why this matters. Every other script adds reflections as POWER (incoherent sum) -
the standard assumption for a swept, comparatively broadband interrogation. That
is only valid if the source coherence length is short compared with the optical
path differences between reflectors. The procured bench spaces gratings ~4 m apart
(round-trip path difference ~8 m in fiber), so the question is whether L_coh from
the VCSEL line is well below that. If it is, power addition is justified and the
coherent-crosstalk term in this script is an over-pessimistic bound, not a real
risk. If L_coh were comparable to the spacing, the field-sum spread below would be
the actual added error.

Coherence length. For a Lorentzian line of FWHM delta_nu the (1/e field) coherence
length is L_coh = c / (pi * delta_nu). We evaluate it for the CW VCSEL line
(common.LASER_LINEWIDTH_GHZ ~ 0.3 GHz) and for a broadened ~GHz line, both in
FIBER (c/n, n = 1.468), and compare to the 4 m grating spacing.

Bench case. Two co-located IDENTICAL gratings at the same wavelength (the worst
case for CDM crosstalk: same code residue, zero delay separation). We recover the
apparent Bragg wavelength when their reflections add
  (a) INCOHERENTLY:  |E1|^2 + |E2|^2                (power sum, paper's assumption)
  (b) COHERENTLY:    |E1 + E2 e^{i*phi}|^2          (field sum, phi swept 0..2pi)
and report the peak-to-peak spread of the recovered wavelength over phi. That
spread is the extra error coherent beating would add.
"""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore'); import common as C
np.random.seed(9)

PM = C.PM_PER_GHZ
F  = C.FBG_FWHM_GHZ
R  = C.FBG_R
C_LIGHT = 2.99792458e8             # m/s
N_FIBER = 1.468                    # silica group index @1550 nm
SPACING_M = 4.0                    # grating spacing on the bench

def coherence_length_m(delta_nu_ghz):
    """L_coh = c/(pi*delta_nu) in FIBER (Lorentzian, 1/e field coherence)."""
    return (C_LIGHT / N_FIBER) / (np.pi * delta_nu_ghz * 1e9)

# ---------------------------------------------------------------------------
# coherence length: CW line vs a broadened line
# ---------------------------------------------------------------------------
dnu_cw   = C.LASER_LINEWIDTH_GHZ          # ~0.30 GHz
dnu_broad = 1.0                            # broadened line ~1 GHz (e.g. dithered)
Lcoh_cw   = coherence_length_m(dnu_cw)
Lcoh_broad = coherence_length_m(dnu_broad)

# ---------------------------------------------------------------------------
# two co-located identical gratings: incoherent (power) vs coherent (field) sum
# ---------------------------------------------------------------------------
# Field reflection of a single grating ~ sqrt(power lineshape). tanh^2 power ->
# |r(nu)| = sqrt(R) * sech((nu-nu_b)/sig). We take both gratings identical and
# co-located (same nu_b, negligible spectral offset) - the hardest crosstalk case.
# Grid is fine (step ~ F/265) so the fitted peak is not quantization-limited.
nu = np.linspace(-6*F, 6*F, 3200)
nu_b = 0.0
amp = np.sqrt(R) * np.sqrt(np.clip(C.fbg_tanh(nu, nu_b, F), 0, None))   # field magnitude

# (a) incoherent power sum
P_incoh = amp**2 + amp**2                    # = 2*R*sech^2 ; peak unchanged
lam_incoh = C.gauss_fit_peak(nu, P_incoh) * PM

# (b) coherent field sum with relative phase phi swept 0..2pi
phis = np.linspace(0.0, 2*np.pi, 73)
lam_coh = np.empty_like(phis)
for i, phi in enumerate(phis):
    E = amp + amp * np.exp(1j*phi)           # co-located, identical field shapes
    P = np.abs(E)**2                          # |E1 + E2 e^{iphi}|^2
    lam_coh[i] = C.gauss_fit_peak(nu, P) * PM

coh_pp   = float(np.ptp(lam_coh))            # peak-to-peak wavelength spread [pm]
coh_rms  = float(np.std(lam_coh))            # RMS over phase [pm]
coh_maxdev = float(np.max(np.abs(lam_coh - lam_incoh)))

# A slightly DETUNED pair moves the apparent peak more than the perfectly
# co-located one (the fringe becomes asymmetric about grating 1). Scan the
# separation ONLY over the UNRESOLVED regime (<= 0.5*FWHM), where the two
# reflections blend into a single apparent Bragg peak and the single-peak fit is
# the right model; wider separations resolve into two peaks (a different problem,
# handled by CDM/despreading), so a single-peak "spread" there is not meaningful.
seps = np.array([0.1, 0.2, 0.3, 0.4, 0.5]) * F       # separations [GHz], unresolved
det_pp_by_sep = []
for sep in seps:
    amp2 = np.sqrt(R) * np.sqrt(np.clip(C.fbg_tanh(nu, sep, F), 0, None))
    lam = np.array([C.gauss_fit_peak(nu, np.abs(amp + amp2*np.exp(1j*phi))**2) * PM
                    - C.gauss_fit_peak(nu, amp**2 + amp2**2) * PM for phi in phis])
    det_pp_by_sep.append(float(np.ptp(lam)))
det_pp_by_sep = np.array(det_pp_by_sep)
det_pp = float(det_pp_by_sep.max())                  # worst-case coherent excursion
det_sep_worst = float(seps[int(np.argmax(det_pp_by_sep))] * PM)   # pm
# keep a representative detuned trace (worst separation) for the figure
amp2 = np.sqrt(R) * np.sqrt(np.clip(C.fbg_tanh(nu, seps[int(np.argmax(det_pp_by_sep))], F), 0, None))
lam_coh_det = np.array([C.gauss_fit_peak(nu, np.abs(amp + amp2*np.exp(1j*phi))**2) * PM
                        for phi in phis])

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(10.8, 3.8))
# (a) a few coherent lineshapes vs the incoherent one
ax[0].plot(nu*PM/1000, P_incoh/P_incoh.max(), 'k', lw=1.8, label='incoherent (power sum)')
for phi, c in [(0.0,'#2980b9'), (np.pi/2,'#e67e22'), (np.pi,'#c0392b')]:
    P = np.abs(amp + amp*np.exp(1j*phi))**2
    ax[0].plot(nu*PM/1000, P/P_incoh.max(), '--', lw=1.0,
               label=f'coherent, phi={phi/np.pi:.1f}$\\pi$')
ax[0].set_xlabel('wavelength offset [nm]'); ax[0].set_ylabel('norm. reflected power')
ax[0].set_title('(a) Co-located identical pair: field sum vs power sum')
ax[0].legend(fontsize=7)
# (b) recovered Bragg wavelength vs relative phase
ax[1].axhline(lam_incoh, color='k', lw=1.4, label='incoherent (paper)')
ax[1].plot(phis/np.pi, lam_coh, 'o-', ms=3, color='#2980b9',
           label=f'coherent, identical (pp {coh_pp:.2f} pm)')
ax[1].plot(phis/np.pi, lam_coh_det, 's-', ms=3, color='#c0392b',
           label=f'coherent, detuned {det_sep_worst:.0f} pm (pp {det_pp:.2f} pm)')
ax[1].set_xlabel('relative phase $\\phi/\\pi$'); ax[1].set_ylabel('recovered $\\lambda_B$ [pm]')
ax[1].set_title('(b) Apparent Bragg wavelength vs beating phase')
ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig('figs/fig_s9_coherent.png', dpi=140)

# ---------------------------------------------------------------------------
# printed numbers + verdict
# ---------------------------------------------------------------------------
print(f"source line (CW, Markowski 2023): delta_nu = {dnu_cw:.2f} GHz")
print(f"  L_coh (CW,   in fiber n={N_FIBER}) = {Lcoh_cw:.3f} m")
print(f"  L_coh (broadened ~{dnu_broad:.0f} GHz)   = {Lcoh_broad:.3f} m")
print(f"  grating spacing on bench          = {SPACING_M:.1f} m  "
      f"(round-trip path diff ~ {2*SPACING_M:.0f} m)")
print(f"  ratio L_coh(CW) / spacing         = {Lcoh_cw/SPACING_M:.3f}")
print("--- two co-located identical gratings (worst-case crosstalk) ---")
print(f"  incoherent (power sum) recovered lambda_B = {lam_incoh:+.3f} pm")
print(f"  coherent (field sum, phi 0..2pi): peak-to-peak spread = {coh_pp:.3f} pm, "
      f"RMS = {coh_rms:.3f} pm")
print(f"  max deviation from incoherent             = {coh_maxdev:.3f} pm")
print(f"  (identical/co-located: symmetric beating leaves the peak nearly fixed)")
print(f"  detuned pair, separation scan (0.1..0.5*FWHM, unresolved regime):")
for sep, pp in zip(seps, det_pp_by_sep):
    print(f"     sep = {sep*PM:6.0f} pm ({sep/F:.2f}*FWHM): peak-to-peak = {pp:.3f} pm")
print(f"  worst-case coherent excursion = {det_pp:.2f} pm at {det_sep_worst:.0f} pm separation  "
      f"<-- the coherent-risk bound")
print("--- VERDICT (coherent vs incoherent regime) ---")
if Lcoh_cw < 0.5*SPACING_M:
    print(f"  L_coh(CW) = {Lcoh_cw:.2f} m is well below the {SPACING_M:.0f} m grating spacing, so "
          f"reflections from different gratings arrive mutually incoherent and add in POWER; "
          f"the paper's incoherent assumption holds and the co-located field-sum spread "
          f"(<= {max(coh_pp,det_pp):.1f} pm) is a conservative worst-case bound, not the expected error.")
    print(f"  The one place beating can bite is two reflectors within ~L_coh of each other "
          f"(e.g. a grating and a nearby connector), not the 4 m-spaced sensor array.")
else:
    print(f"  L_coh(CW) = {Lcoh_cw:.2f} m is comparable to the {SPACING_M:.0f} m spacing: the bench "
          f"may be partially coherent and the field-sum spread up to ~{max(coh_pp,det_pp):.1f} pm "
          f"should be treated as a real error term.")
print("saved figs/fig_s9_coherent.png")
