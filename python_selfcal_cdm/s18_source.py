"""s18_source.py - where the source hurts: chirp of the code, and what it does
on the edge of a grating.

This is the explanatory figure for the term that gives the paper its name, and it
was missing from the earlier draft. Three panels.

  (a) What "chirp of the code" means. Modulating the laser current does not only
      switch the power; it moves the optical frequency. Two components appear:
      an adiabatic shift that follows the instantaneous power, so the laser sits
      at a different frequency during an "on" chip than during an "off" chip, and
      a transient overshoot at every edge, following the derivative of the power.
      The laser therefore does not probe a single optical frequency while the
      code runs; it probes a distribution. The inset is that distribution, and it
      is exactly the kernel p(delta) that the model integrates over. Everything
      is drawn in units of the excursion Delta, because the excursion of the
      BW10 VCSEL under code modulation has not been measured.

  (b) Why the distribution matters. On the flank of a grating the reflectance
      changes steeply with optical frequency, so averaging over the chirp
      distribution is not the same as sampling at its mean: a skewed kernel on an
      asymmetric flank returns a biased value. Summed across the sweep, the
      recovered line is displaced. That is the FM-to-AM error.

  (c) How the references calibrate. At one operating point (0.3 FWHM of
      chirp, plus a drift of the lambda(V) table of +20 pm and 3 percent
      gain) the axis error across the band, what each reference reads at
      its own band position, and the polynomial through the readings: a
      constant from one reference, a line from two, a parabola from three.
      A single reference sits at the sweep centre. The panel spans the
      whole 10-nm sweep: three references in the edge bands and the centre
      band serve every band, because the axis error is smooth over the sweep.

  (d) How far the references get. The chirp offset varies smoothly across the
      band, so temperature-stabilised references sample it and a low-order fit
      removes it. One reference removes the constant, two the slope, three the
      curvature. What survives scales with the excursion and with the spread of
      grating lineshapes.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS
FS.apply()


PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
BAND_HALF = 625.0                    # +-5 nm, the whole VCSEL sweep. The axis
                                     # error is smooth over the sweep, so the
                                     # references are spread over it, not over
                                     # one 400-pm band. The residual does not
                                     # depend on this width (normalised).

# ---------------------------------------------------------------------------
# (a) chirp waveform during code modulation, in units of the excursion
# ---------------------------------------------------------------------------
OSR = 64
chips = np.array([1, 0, 1, 1, 0, 1, 0, 0], float)
drive = np.repeat(chips, OSR)
t = np.arange(len(drive)) / float(OSR)          # time in chip periods

# adiabatic part: the frequency follows the instantaneous power with a short lag
tau_ad = 0.02
ad = np.zeros_like(drive)
for i in range(1, len(drive)):
    ad[i] = ad[i - 1] + (drive[i] - ad[i - 1]) / (tau_ad * OSR)
# transient part: follows the derivative of the power and damps out over a
# fraction of a chip, so it rides on top of the adiabatic level
d = np.gradient(drive) * OSR
tr = np.zeros_like(drive)
decay = np.exp(-1.0 / (0.10 * OSR))
for i in range(1, len(drive)):
    tr[i] = tr[i - 1] * decay + d[i] / OSR
tr = tr / max(np.abs(tr).max(), 1e-9)
nu_raw = ad + 0.35 * tr
plateau = float(np.median(nu_raw[drive > 0.5]))
nu_inst = (nu_raw - nu_raw.min()) / (nu_raw.max() - nu_raw.min())   # delta_nu/Delta
plateau = (plateau - nu_raw.min()) / (nu_raw.max() - nu_raw.min())

# the kernel: distribution of instantaneous frequency, weighted by optical power
hist, edges = np.histogram(nu_inst, bins=36, weights=drive, density=True)
centres = 0.5 * (edges[1:] + edges[:-1])

# ---------------------------------------------------------------------------
# (b) FM-to-AM on the grating flank
# ---------------------------------------------------------------------------
DELTA_DEMO = 0.30 * F                # excursion used for the illustration only
ASYM = 0.30
g = np.linspace(-2.6 * F, 2.6 * F, 1200)
true_line = C.fbg_tanh(g, 0.0, F, n_side=ASYM)
chirped = C.fmam_readout(g, 0.0, F, DELTA_DEMO, mean_off=0.30 * DELTA_DEMO,
                         skew=1.2, shape='tanh', asym=ASYM)
chirped = chirped / chirped.max()
p_true = C.gauss_fit_peak(g, true_line) * PM
p_chirp = C.gauss_fit_peak(g, chirped) * PM

# ---------------------------------------------------------------------------
# (c) residual against excursion and reference count
# ---------------------------------------------------------------------------
def chirp_residual(ratio, nref, nsen=8, seed=3):
    rng = np.random.default_rng(seed)
    delta = ratio * F
    ref_nu = (np.linspace(-0.9 * BAND_HALF, 0.9 * BAND_HALF, nref) if nref > 1
              else np.array([0.0]))
    sen_nu = np.sort(rng.uniform(-BAND_HALF, BAND_HALF, nsen))
    sen_as = rng.uniform(-0.30, 0.30, nsen)
    ref_as = rng.uniform(-0.05, 0.05, max(nref, 1))

    def off(nb):
        x = nb / BAND_HALF
        return delta * (0.20 + 0.20 * x + 0.15 * x ** 2 + 0.10 * np.sin(2.5 * x))

    def shift(nb, asym):
        gg = np.linspace(nb - 5 * F, nb + 5 * F, 801)
        a = C.fmam_readout(gg, nb, F, delta, mean_off=off(nb), skew=1.2,
                           shape='tanh', asym=asym)
        b = C.fbg_tanh(gg, nb, F, n_side=asym)
        return (C.gauss_fit_peak(gg, a) - C.gauss_fit_peak(gg, b)) * PM

    se = np.array([shift(sen_nu[i], sen_as[i]) for i in range(nsen)])
    if nref == 0:
        return float(np.sqrt(np.mean(se ** 2)))
    re = np.array([shift(ref_nu[j], ref_as[j]) for j in range(nref)])
    p = np.polyfit(ref_nu, re, min(nref - 1, 2))
    return float(np.sqrt(np.mean((se - np.polyval(p, sen_nu)) ** 2)))


ratios = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1.0])

# ---------------------------------------------------------------------------
# (c) the calibration itself at one operating point
# ---------------------------------------------------------------------------
CAL_RATIO, DRIFT_OFF, DRIFT_GAIN = 0.30, 20.0, 0.003  # FWHM, pm, relative
# drift of the lambda(V) table since the last uncoded calibration: 0.2 K of
# laser temperature at 0.1 nm/K is 20 pm of offset, and 1e-4/K of spring
# stiffness over 30 K is a 0.3 percent gain error, 15 pm at the sweep edges


def calibration(ratio=CAL_RATIO, nsen=8, seed=3):
    rng = np.random.default_rng(seed)
    delta = ratio * F
    sen_nu = np.sort(rng.uniform(-BAND_HALF, BAND_HALF, nsen))
    sen_as = rng.uniform(-0.30, 0.30, nsen)

    def off(nb):
        x = nb / BAND_HALF
        return delta * (0.20 + 0.20 * x + 0.15 * x ** 2 + 0.10 * np.sin(2.5 * x))

    def chirp_shift(nb, asym):
        gg = np.linspace(nb - 5 * F, nb + 5 * F, 801)
        a = C.fmam_readout(gg, nb, F, delta, mean_off=off(nb), skew=1.2,
                           shape='tanh', asym=asym)
        b = C.fbg_tanh(gg, nb, F, n_side=asym)
        return (C.gauss_fit_peak(gg, a) - C.gauss_fit_peak(gg, b)) * PM

    def drift(nb):
        return DRIFT_OFF + DRIFT_GAIN * nb * PM

    def err(nb, asym):
        return drift(nb) + chirp_shift(nb, asym)

    grid = np.linspace(-BAND_HALF, BAND_HALF, 121)
    out = dict(grid=grid, drift=drift(grid),
               smooth=np.array([err(x, 0.0) for x in grid]),
               sen_nu=sen_nu,
               sen_err=np.array([err(sen_nu[i], sen_as[i]) for i in range(nsen)]),
               refs={})
    for nref in (1, 2, 3):
        ref_nu = (np.linspace(-0.9 * BAND_HALF, 0.9 * BAND_HALF, nref) if nref > 1
                  else np.array([0.0]))
        ref_as = rng.uniform(-0.05, 0.05, nref)
        ref_rd = np.array([err(ref_nu[j], ref_as[j]) for j in range(nref)])
        coef = np.polyfit(ref_nu, ref_rd, nref - 1)
        res = out['sen_err'] - np.polyval(coef, sen_nu)
        out['refs'][nref] = dict(nu=ref_nu, rd=ref_rd, fit=np.polyval(coef, grid),
                                 res=res, rms=float(np.sqrt(np.mean(res ** 2))))
    return out


cal = calibration()
curves = {n: np.array([chirp_residual(r, n) for r in ratios]) for n in (0, 1, 2, 3)}

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
# The paper places this at full text width (7.16 in). Draw at the final
# physical width so annotations remain legible after inclusion.
fig, ax = plt.subplots(1, 4, figsize=(7.16, 1.95))

# --- (a) -------------------------------------------------------------------
ax[0].fill_between(t, -0.4, 1.6, where=drive > 0.5, step='post',
                   color='0.88', lw=0)
ax[0].plot(t, nu_inst, color='#D55E00', lw=1.4)
ax[0].set_ylim(-0.18, 1.45)
ax[0].set_xlim(0, 8)
ax[0].set_xlabel('time [chip periods]')
ax[0].set_ylabel('code-induced shift  ' + r'$\delta\lambda/\Delta\lambda_{\mathrm{ch}}$')
FS.letter(ax[0], 'a')
ax[0].axhline(plateau, color='#D55E00', ls=':', lw=0.8)
ax[0].annotate('', xy=(1.62, 0.0), xytext=(1.62, 1.0),
               arrowprops=dict(arrowstyle='<->', lw=0.9, color='#0072B2'))
ax[0].text(1.72, 0.5, r'$\Delta\lambda_{\mathrm{ch}}$', color='#0072B2', fontsize=8.0, ha='left',
           va='center')
ax[0].grid(False, alpha=0.2)

# --- (b) -------------------------------------------------------------------
ax[1].plot(g * PM / 1000.0, true_line, color='#0072B2', lw=1.4, label='$R_k(\\lambda)$, no chirp')
ax[1].plot(g * PM / 1000.0, chirped, color='#D55E00', lw=1.4,
           label='$S_k^{\\mathrm{ch}}$, chirped source')
nu_op = -0.62 * F
kern = np.exp(-0.5 * ((g - nu_op) / (0.5 * DELTA_DEMO)) ** 2)
ax[1].fill_between(g * PM / 1000.0, 0, 0.34 * kern, color='#CC79A7', alpha=0.32,
                   lw=0)
ax[1].annotate(r'$p(\xi)$', xy=(nu_op * PM / 1000.0, 0.30),
               xytext=(-0.52, 0.62), fontsize=7.4, color='#CC79A7',
               ha='center', va='center',
               arrowprops=dict(arrowstyle='-', color='#CC79A7', lw=0.6))
ax[1].axvline(p_true / 1000.0, color='#0072B2', ls=':', lw=0.9)
ax[1].axvline(p_chirp / 1000.0, color='#D55E00', ls=':', lw=0.9)
# 67 pm to piec procent szerokosci panelu, wiec groty ida na zewnatrz
# linii wymiarowych i pokazuja do srodka, a liczba przenosi sie obok
FS.dim_gap(ax[1], p_true / 1000.0, p_chirp / 1000.0, 1.10,
           '%.0f pm' % abs(p_chirp - p_true), color='0.2',
           tail=0.13, side='left')
ax[1].set_ylim(0, 1.62)
ax[1].set_xlim(-0.62, 0.62)
ax[1].set_xlabel('wavelength offset [nm]')
ax[1].set_ylabel('normalized reflectance')
FS.letter(ax[1], 'b')
ax[1].legend(fontsize=5.8, loc='upper right',
             ncol=1, frameon=True, handlelength=1.5, columnspacing=0.8,
             labelspacing=0.18, borderaxespad=0.25)

# --- (c) -------------------------------------------------------------------
GREY = '0.45'
xpm = cal['grid'] * PM / 1000.0
ax[2].axhline(0, color='0.6', lw=0.6)
chirp_part = cal['smooth'] - cal['drift']
ax[2].plot(xpm, cal['drift'], color='0.25', lw=0.8, ls=':')
ax[2].plot(xpm, cal['smooth'], color='0.25', lw=1.3)
ax[2].text(4.8, cal['drift'][-1] + 2.5, 'drift of the $\\lambda(V)$ table', fontsize=5.2,
           color='0.25', ha='right', va='bottom')
i60 = int(np.argmin(np.abs(xpm + 1.5)))
ax[2].text(-1.5, cal['smooth'][i60] + 3.5, 'axis error: drift + chirp', fontsize=5.2,
           color='0.25', ha='center', va='bottom')
ax[2].plot(cal['sen_nu'] * PM / 1000.0, cal['sen_err'], 'o', color=GREY, ms=3.0, mfc='white',
           mew=0.8, label='sensors')
cstyle = {1: ('#E69F00', 's', '1 ref'), 2: ('#0072B2', '^', '2 refs'), 3: ('#009E73', 'd', '3 refs')}
for n in (1, 2, 3):
    col, mk, lab = cstyle[n]
    r = cal['refs'][n]
    ax[2].plot(xpm, r['fit'], color=col, lw=0.9, ls='--')
    ax[2].plot(r['nu'] * PM / 1000.0, r['rd'], mk, color=col, ms=4.2, label=lab + ', fit')
    ax[2].plot(cal['sen_nu'] * PM / 1000.0, r['res'], mk, color=col, ms=2.4, mew=0)
ax[2].plot([], [], 'o', color='0.3', ms=2.4, mew=0, label='after fit')
ax[2].text(4.8, 4.0, 'after the fit', fontsize=5.2, color='0.25', ha='right', va='bottom')
ax[2].set_xlim(-5, 5)
lo = min(cal['smooth'].min(), cal['sen_err'].min())
hi = max(cal['smooth'].max(), cal['drift'].max())
ax[2].set_ylim(lo - 0.55 * (hi - lo), hi + 14)
ax[2].set_xlabel('sweep position [nm]')
ax[2].set_ylabel('reported $-$ true [pm]')
FS.letter(ax[2], 'c')
ax[2].legend(fontsize=5.0, loc='lower left', ncol=2, frameon=True, handlelength=1.4,
             columnspacing=0.6, labelspacing=0.15, borderaxespad=0.25)
ax[2].grid(False, alpha=0.2)

# --- (d) -------------------------------------------------------------------
styles = {0: ('o-', '#D55E00', 'no ref'), 1: ('s-', '#E69F00', '1 ref'),
          2: ('^-', '#0072B2', '2 refs'), 3: ('d-', '#009E73', '3 refs')}
for n in (0, 1, 2, 3):
    mk, col, lab = styles[n]
    ax[3].plot(ratios, curves[n], mk, color=col, lw=1.2, ms=4, label=lab)
ax[3].axhline(10.0, color='0.3', ls='--', lw=0.8, label='10 pm target')
ax[3].set_yscale('log')
ax[3].set_xlabel(r'chirp span  $\Delta\lambda_{\mathrm{ch}}/\mathrm{FWHM}$')
ax[3].set_ylabel('residual $\\delta\\lambda_k$ [pm]')
FS.letter(ax[3], 'd')
ax[3].set_ylim(0.5, 6000)
ax[3].legend(fontsize=5.2, loc='upper left',
             ncol=2, frameon=True, handlelength=1.2, columnspacing=0.5,
             labelspacing=0.15)
ax[3].grid(False, which='both', alpha=0.25)

fig.subplots_adjust(left=0.065, right=0.99, top=0.87, bottom=0.17,
                    wspace=0.40)
fig.savefig('figs/fig_s18_source.png', dpi=150, bbox_inches='tight')
fig.savefig('figs/fig_s18_source.pdf', bbox_inches='tight')

print('illustration only: excursion %.2f x FWHM, asymmetry %.2f, apparent shift %.1f pm'
      % (DELTA_DEMO / F, ASYM, p_chirp - p_true))
print('--- residual [pm] vs excursion and reference count ---')
print('  Delta/FWHM ' + ''.join('%9.2f' % r for r in ratios))
for n in (0, 1, 2, 3):
    print('  %d ref%s     ' % (n, ' ' if n == 1 else 's') +
          ''.join('%9.2f' % v for v in curves[n]))
print('(c) at %.2f FWHM, drift %+.0f pm and %.1f%% gain: RMS after 1/2/3 refs = %s pm'
      % (CAL_RATIO, DRIFT_OFF, 100 * DRIFT_GAIN,
         [round(cal['refs'][n]['rms'], 2) for n in (1, 2, 3)]))
print('saved figs/fig_s18_source.png')
