"""s54_refcal.py - how the reference gratings calibrate the wavelength axis.

Companion to s18 (c). That panel gives only the residual against the chirp
span. This one shows the mechanism at one operating point:

  (a) The axis error across the band: what the interrogator reads minus the
      true Bragg wavelength. It has two parts, the drift of the lambda(V)
      table since the last uncoded calibration (an offset and a gain) and the
      chirp of the code (the s18 kernel, whose mean offset varies smoothly
      across the band). Reference gratings at a known, stabilised wavelength
      read that error at their own band positions. A polynomial through the
      readings is the correction: a constant from one reference, a line from
      two, a parabola from three.
  (b) What is left on the sensors after subtracting each polynomial.

Same chirp model and same sensor set as s18, chirp span 0.3 FWHM, drift
+20 pm offset and 3 percent gain. The sensors carry random line asymmetries,
so their chirp shift differs a little from the smooth curve, which is the
part no fit can remove.

Output: figs/fig_s54_refcal.pdf/png, column width.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS
FS.apply()

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
BAND_HALF = 25.0                     # +-200 pm, the band W of Table III
RATIO = 0.30                         # chirp span / FWHM
DRIFT_OFF = 20.0                     # pm, table offset since calibration
DRIFT_GAIN = 0.03                    # relative gain error of the table
NSEN = 8
SEED = 3

BLUE, VERM, ORAN, GREE = '#0072B2', '#D55E00', '#E69F00', '#009E73'
GREY = '0.45'

rng = np.random.default_rng(SEED)
delta = RATIO * F / 2.0
sen_nu = np.sort(rng.uniform(-BAND_HALF, BAND_HALF, NSEN))
sen_as = rng.uniform(-0.30, 0.30, NSEN)


def off(nb):
    """mean chirp offset of the kernel across the band, as in s18"""
    x = nb / BAND_HALF
    return delta * (0.20 + 0.20 * x + 0.15 * x ** 2 + 0.10 * np.sin(2.5 * x))


def chirp_shift(nb, asym):
    """fitted-centre shift of one grating under the chirp kernel, pm"""
    gg = np.linspace(nb - 5 * F, nb + 5 * F, 801)
    a = C.fmam_readout(gg, nb, F, delta, mean_off=off(nb), skew=1.2,
                       shape='tanh', asym=asym)
    b = C.fbg_tanh(gg, nb, F, n_side=asym)
    return (C.gauss_fit_peak(gg, a) - C.gauss_fit_peak(gg, b)) * PM


def drift(nb):
    """drift of the lambda(V) table: offset plus gain, pm"""
    return DRIFT_OFF + DRIFT_GAIN * nb * PM


def axis_error(nb, asym):
    return drift(nb) + chirp_shift(nb, asym)


# smooth curve: symmetric lines, so only the smooth part
grid = np.linspace(-BAND_HALF, BAND_HALF, 121)
smooth = np.array([axis_error(g, 0.0) for g in grid])
sen_err = np.array([axis_error(sen_nu[i], sen_as[i]) for i in range(NSEN)])

refs = {}
for nref in (1, 2, 3):
    ref_nu = np.linspace(-0.9 * BAND_HALF, 0.9 * BAND_HALF, nref) if nref > 1 \
        else np.array([0.0])
    ref_as = rng.uniform(-0.05, 0.05, nref)
    ref_rd = np.array([axis_error(ref_nu[j], ref_as[j]) for j in range(nref)])
    coef = np.polyfit(ref_nu, ref_rd, nref - 1)
    refs[nref] = dict(nu=ref_nu, rd=ref_rd, coef=coef,
                      fit=np.polyval(coef, grid),
                      res=sen_err - np.polyval(coef, sen_nu))

rms = {0: float(np.sqrt(np.mean(sen_err ** 2)))}
for n in (1, 2, 3):
    rms[n] = float(np.sqrt(np.mean(refs[n]['res'] ** 2)))

# ---------------------------------------------------------------------------
fig, ax = plt.subplots(2, 1, figsize=(3.45, 4.4), sharex=True)
fig.subplots_adjust(left=0.16, right=0.98, top=0.95, bottom=0.095, hspace=0.30)
xpm = grid * PM
chirp_only = smooth - drift(grid)

a = ax[0]
a.axhline(0, color='0.6', lw=0.6)
a.plot(xpm, drift(grid), color='0.25', lw=0.8, ls=':', label='drift of the $\\lambda(V)$ table')
a.plot(xpm, chirp_only, color='0.25', lw=0.8, ls='-.', label='chirp of the code')
a.plot(xpm, smooth, color='0.25', lw=1.4, label='axis error, their sum')
a.plot(sen_nu * PM, sen_err, 'o', color=GREY, ms=3.5, mfc='white', mew=0.9,
       label='sensors, as read')
styles = {1: (ORAN, 's', 'one reference'), 2: (BLUE, '^', 'two references'),
          3: (GREE, 'D', 'three references')}
for n in (1, 2, 3):
    col, mk, lab = styles[n]
    a.plot(xpm, refs[n]['fit'], color=col, lw=1.0, ls='--')
    a.plot(refs[n]['nu'] * PM, refs[n]['rd'], mk, color=col, ms=5.0,
           label=lab + ', fit')
a.set_ylabel('reported $-$ true $\\lambda_B$ [pm]')
a.legend(fontsize=5.4, loc='lower left', ncol=2, handlelength=1.8,
         columnspacing=0.8)
lo = min(smooth.min(), sen_err.min())
hi = max(smooth.max(), sen_err.max(), drift(grid).max())
a.set_ylim(lo - 0.55 * (hi - lo) - 5, hi + 8)
FS.letter(a, 'a')

b = ax[1]
b.axhspan(-10, 10, color='0.92', lw=0, zorder=0)
b.axhline(0, color='0.6', lw=0.6)
for n in (1, 2, 3):
    col, mk, lab = styles[n]
    b.plot(sen_nu * PM, refs[n]['res'], mk, color=col, ms=4.5, ls='-', lw=0.6,
           label='%s, %.1f pm RMS' % (lab, rms[n]))
b.plot([], [], ' ', label='no correction, %.0f pm RMS' % rms[0])
b.text(198, 10.8, '$\\pm$10 pm', fontsize=5.6, color='0.4', va='bottom', ha='right')
b.set_xlabel('band position [pm]')
b.set_ylabel('sensor error after correction [pm]')
b.set_xlim(-200, 200)
b.set_ylim(-22, 22)
b.legend(fontsize=5.6, loc='lower left', ncol=1, handlelength=1.8)
FS.letter(b, 'b')

for ext in ('pdf', 'png'):
    fig.savefig('figs/fig_s54_refcal.' + ext, dpi=300)

print('chirp span %.2f FWHM = %.0f pm, drift %+.0f pm and %.0f%% gain'
      % (RATIO, delta * PM, DRIFT_OFF, 100 * DRIFT_GAIN))
print('sensor positions [pm]:', np.round(sen_nu * PM, 0))
print('axis error on sensors [pm]:', np.round(sen_err, 1))
for n in (1, 2, 3):
    print('%d ref: readings %s pm, coef %s, residual %s pm'
          % (n, np.round(refs[n]['rd'], 1), np.round(refs[n]['coef'], 3),
             np.round(refs[n]['res'], 1)))
print('RMS [pm]: none %.2f, 1 ref %.2f, 2 refs %.2f, 3 refs %.2f'
      % (rms[0], rms[1], rms[2], rms[3]))
print('saved figs/fig_s54_refcal.pdf')
