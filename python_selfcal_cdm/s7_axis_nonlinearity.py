"""s7_axis_nonlinearity.py - the nonlinear voltage->wavelength tuning of the
HCG-VCSEL, and how co-coded reference gratings straighten that axis.

The interrogator ramps the tuning voltage v over [0,1] uniformly across M slots.
The TRUE optical frequency is nonlinear in v (HCG-VCSEL ~ quadratic, BW10-1550):
    nu_true(v) = common.vcsel_sweep(v, span, quad=0.30).
A naive interrogator instead ASSUMES a straight line over the same span,
    nu_lin(v) = v * span.
A grating with true Bragg frequency nu_B reflects at the voltage v* where
nu_true(v*) = nu_B; the naive readout reports nu_lin(v*). The axis error,
nu_lin(v*) - nu_B, is the sweep nonlinearity biasing the apparent Bragg position.

Reference gratings at temperature-set, KNOWN frequencies give anchor points
(v*, nu_B). Fitting the v->nu map through them calibrates the axis:
  - 2 references -> LINEAR fit (cannot reproduce the HCG curvature -> residual);
  - 3 references -> QUADRATIC fit (captures the dominant curvature; a small
    higher-order residual remains);
  - an MZI k-clock supplies the dense/true map and clears even that residual.
On top of the axis sits the source-chirp FM-AM error (as in s3/s6): an apparent
shift per grating from common.fmam_readout at the measured chirp. The MZI is a
source-frequency ruler and is blind to that detection-side effect, so FM-AM
still needs the references. Total per-sensor error = axis error + FM-AM error.

Teaching point: the minimal 3-grating bench (2 ref + 1 sensor) only fits a
LINEAR axis and leaves residual HCG curvature; 3 references (from the 15
procured) or the MZI remove it, while FM-AM removal always needs the references.
"""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore'); import common as C

PM   = C.PM_PER_GHZ
F    = C.FBG_FWHM_GHZ          # ~31.25 GHz (= 250 pm)
SPAN = 600.0                   # general-case sweep band [GHz] (~4.8 nm)
QUAD = 0.30                    # HCG-VCSEL quadratic tuning nonlinearity (dominant)
CUBE = 0.10                    # small higher-order (cubic) term: a quadratic fit
                               # cannot fully reproduce it, only the dense MZI can
CHIRP = C.MEAS_CHIRP_GHZ       # measured source chirp ~10.6 GHz (FM-AM driver)
M    = 2048                    # voltage slots in the ramp (dense grid)

# --- the sweep on a fine voltage grid: true (quadratic) vs naive (linear) ----
v_grid   = np.linspace(0.0, 1.0, M)
nu_true  = C.vcsel_sweep(v_grid, SPAN, quad=QUAD, cube=CUBE)   # nonlinear ground truth
nu_lin   = v_grid * SPAN                             # naive straight-line map
# nu_true is monotone in v -> invert it to get the voltage that hits any nu_B
def v_star(nu_b):
    return float(np.interp(nu_b, nu_true, v_grid))
# what the naive interrogator REPORTS for a grating whose true freq is nu_b:
def nu_naive(nu_b):
    return float(np.interp(v_star(nu_b), v_grid, nu_lin))

# --- FM-AM apparent shift of one grating (same recipe as s3/s6) --------------
def fmam_err(nu_b, asym, skew=1.2):
    """Apparent Bragg shift [GHz] from source-chirp FM->AM at the measured chirp."""
    g = np.linspace(nu_b - 6*F, nu_b + 6*F, 1400)
    spec = C.fmam_readout(g, nu_b, F, CHIRP, mean_off=0.30*CHIRP,
                          skew=skew, shape='tanh', asym=asym)
    return C.gauss_fit_peak(g, spec) - nu_b

# --- sensor / reference layout ------------------------------------------------
NSEN = 12                                    # sensor gratings spanning the band
NREF = 3                                      # up to 3 co-coded references
NR   = 30                                     # random realizations
# references placed across the band (low / mid / high) at known set-points:
ref_nu = np.array([0.10, 0.50, 0.90]) * SPAN  # true frequencies (temperature-set)
mae = lambda e: np.mean(np.abs(e)) * PM       # mean abs error in pm

# Calibrate the v->nu axis from K references: fit nu_B vs measured voltage v*,
# evaluate at the sensor voltages. order = K-1 (2 ref -> linear, 3 ref -> quad).
def axis_calibrate(sen_vstar, ref_vstar, ref_true, order):
    p = np.polyfit(ref_vstar, ref_true, order)
    return np.polyval(p, sen_vstar)           # calibrated frequency at each sensor

acc = np.zeros(5)
last = None
rng = np.random.default_rng(7)
for it in range(NR):
    # sensors: random positions across the band, varied lineshape asymmetry
    sen_nu   = np.sort(rng.uniform(0.05, 0.95, NSEN)) * SPAN
    sen_asym = rng.uniform(-0.30, 0.30, NSEN)         # sensors vary
    ref_asym = rng.uniform(-0.05, 0.05, NREF)         # references well characterized

    # measured reflection voltages (true), and the naive reported frequencies
    sen_vstar = np.array([v_star(n) for n in sen_nu])
    ref_vstar = np.array([v_star(n) for n in ref_nu])
    sen_naive = np.array([nu_naive(n) for n in sen_nu])   # naive freq readout

    # axis-nonlinearity error of each sensor (naive report minus truth)
    ax_sen = sen_naive - sen_nu

    # FM-AM apparent shift per grating
    fm_sen = np.array([fmam_err(sen_nu[i], sen_asym[i]) for i in range(NSEN)])
    fm_ref = np.array([fmam_err(ref_nu[j], ref_asym[j]) for j in range(NREF)])

    # (1) no calibration: both errors raw
    e_nocal = ax_sen + fm_sen

    # (2) 2 references -> LINEAR axis fit (leaves residual HCG curvature),
    #     and the same 2 refs subtract FM-AM in common mode (linear)
    nu_cal2 = axis_calibrate(sen_vstar, ref_vstar[:2], ref_nu[:2], 1)
    ax_res2 = nu_cal2 - sen_nu                          # residual axis error
    p_fm2   = np.polyfit(ref_nu[:2], fm_ref[:2], 1)
    fm_res2 = fm_sen - np.polyval(p_fm2, sen_nu)
    e_2ref  = ax_res2 + fm_res2

    # (3) 3 references -> QUADRATIC axis fit (captures the dominant HCG quadratic;
    #     a small higher-order axis residual is left),
    #     and a quadratic FM-AM common-mode subtraction
    nu_cal3 = axis_calibrate(sen_vstar, ref_vstar, ref_nu, 2)
    ax_res3 = nu_cal3 - sen_nu
    p_fm3   = np.polyfit(ref_nu, fm_ref, 2)
    fm_res3 = fm_sen - np.polyval(p_fm3, sen_nu)
    e_3ref  = ax_res3 + fm_res3

    # (4) MZI only: dense k-clock fixes the AXIS (true map recovered), but it is
    #     blind to the grating-edge FM-AM -> that error remains uncorrected
    e_mzi   = fm_sen

    # (5) 3 references + MZI: MZI removes the axis exactly (incl. higher order),
    #     references remove FM-AM -> only the FM-AM residual is left
    e_refmzi = fm_res3

    acc += np.array([mae(e_nocal), mae(e_2ref), mae(e_3ref), mae(e_mzi), mae(e_refmzi)])
    if it == 0:
        last = (sen_nu, e_nocal, e_2ref, e_3ref, e_refmzi)

vals = acc / NR
labels_full = ['no calibration', '2 ref (linear axis)', '3 ref (quadratic axis)',
               'MZI only (FM-AM left)', '3 ref + MZI']

print(f"band span = {SPAN:.0f} GHz (~{SPAN*PM/1000:.2f} nm), FBG FWHM = {C.FBG_FWHM_PM:.0f} pm "
      f"(= {F:.1f} GHz), HCG quad = {QUAD}, measured chirp = {CHIRP:.1f} GHz")
print(f"MAE over {NSEN} sensors, mean of {NR} realizations:")
for l, v in zip(labels_full, vals):
    print(f"  {l:24s}: {v:6.1f} pm")

# ---------------------------------------------------------------------------
# Figure: (a) apparent error vs band position, (b) MAE bar chart per strategy
# ---------------------------------------------------------------------------
sn, e0, e2, e3, e5 = last
fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.0))

ax[0].axhline(0, c='0.85', lw=0.7)
ax[0].plot(sn*PM/1000, e0*PM, 'o-', ms=4, color='#c0392b', label='no calibration')
ax[0].plot(sn*PM/1000, e2*PM, 's-', ms=4, color='#e67e22', label='2 ref (linear axis)')
ax[0].plot(sn*PM/1000, e3*PM, 'D-', ms=4, color='#2980b9', label='3 ref (quadratic axis)')
ax[0].plot(sn*PM/1000, e5*PM, '^-', ms=4, color='#27ae60', label='3 ref + MZI')
for rn in ref_nu:
    ax[0].axvline(rn*PM/1000, color='0.75', ls=':', lw=0.8)
ax[0].plot(ref_nu*PM/1000, np.zeros(NREF), 'k^', ms=9, label='references (known T)')
ax[0].set_xlabel('position in band [nm]')
ax[0].set_ylabel(r'apparent $\lambda_B$ error [pm]')
ax[0].set_title('(a) Apparent error vs band position (one realization)')
ax[0].legend(fontsize=7, ncol=2)

bar_labels = ['no\ncalibration', '2 ref\n(linear)', '3 ref\n(quadratic)',
              'MZI only\n(FM-AM left)', '3 ref\n+ MZI']
colors = ['#c0392b', '#e67e22', '#2980b9', '#8e44ad', '#27ae60']
ax[1].bar(range(5), vals, color=colors)
ax[1].set_xticks(range(5)); ax[1].set_xticklabels(bar_labels, fontsize=7.5)
ax[1].set_ylabel(f'MAE [pm] (mean of {NR} runs)')
ax[1].set_title('(b) MAE per calibration strategy')
for i, v in enumerate(vals):
    ax[1].text(i, v + max(vals)*0.02, f'{v:.1f}', ha='center', fontsize=8.5)

plt.tight_layout()
plt.savefig('figs/fig_s7_axis.png', dpi=140)
print("saved figs/fig_s7_axis.png")
