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

  (c) How far the references get. The chirp offset varies smoothly across the
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
BAND_HALF = 28.0                     # +-225 pm window of the bench

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
    ref_nu = np.linspace(-0.9 * BAND_HALF, 0.9 * BAND_HALF, max(nref, 1))
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


ratios = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
curves = {n: np.array([chirp_residual(r, n) for r in ratios]) for n in (0, 1, 2, 3)}

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(7.1, 2.35))

# --- (a) -------------------------------------------------------------------
ax[0].fill_between(t, -0.4, 1.6, where=drive > 0.5, step='post',
                   color='0.88', lw=0)
ax[0].plot(t, nu_inst, color='#c0392b', lw=1.4)
ax[0].set_ylim(-0.18, 1.45)
ax[0].set_xlim(0, 8)
ax[0].set_xlabel('time [chip periods]')
ax[0].set_ylabel('optical frequency shift  ' + r'$\delta\nu/\Delta$')
ax[0].set_title('(a) The code chirps the laser', fontsize=9)
ax[0].text(0.15, -0.11, 'shaded: chip on', fontsize=6.6, color='0.45')
ax[0].axhline(plateau, color='#c0392b', ls=':', lw=0.8)
i_ov = 2 * OSR + int(np.argmax(nu_inst[2 * OSR:3 * OSR]))
ax[0].annotate('transient overshoot', xy=(t[i_ov], nu_inst[i_ov]),
               xytext=(0.22, 1.31), fontsize=6.8, color='#c0392b',
               arrowprops=dict(arrowstyle='->', lw=0.7, color='#c0392b'))
ax[0].annotate('adiabatic level', xy=(3.55, plateau),
               xytext=(2.55, 0.30), fontsize=6.8, color='#c0392b',
               arrowprops=dict(arrowstyle='->', lw=0.7, color='#c0392b'))
ax[0].annotate('', xy=(1.62, 0.0), xytext=(1.62, 1.0),
               arrowprops=dict(arrowstyle='<->', lw=0.9, color='#1f4e79'))
ax[0].text(1.72, 0.5, r'$\Delta$', color='#1f4e79', fontsize=9.5, ha='left',
           va='center')
ax[0].grid(True, alpha=0.2)

axk = ax[0].inset_axes([0.615, 0.50, 0.365, 0.37])
axk.set_facecolor('white')
axk.patch.set_alpha(1.0)
axk.set_zorder(10)
for sp in axk.spines.values():
    sp.set_visible(True); sp.set_linewidth(0.7); sp.set_color('0.45')
axk.fill_between(centres, hist, color='#c0392b', alpha=0.35, lw=0)
axk.plot(centres, hist, color='#c0392b', lw=0.9)
axk.set_title('kernel ' + r'$p(\delta)$', fontsize=6.4, pad=2)
axk.tick_params(labelsize=5.5, pad=1)
axk.set_yticks([])
axk.set_xlim(-0.05, 1.12)
# leader: from the trace into the boxed inset
ax[0].annotate('', xy=(4.66, 1.02), xytext=(3.92, plateau + 0.02),
               arrowprops=dict(arrowstyle='->', lw=0.8, color='0.45',
                               connectionstyle='arc3,rad=0.25'),
               zorder=11, annotation_clip=False)

# --- (b) -------------------------------------------------------------------
ax[1].plot(g * PM / 1000.0, true_line, color='#2980b9', lw=1.4, label='no chirp')
ax[1].plot(g * PM / 1000.0, chirped, color='#c0392b', lw=1.4,
           label='chirped source')
nu_op = -0.62 * F
kern = np.exp(-0.5 * ((g - nu_op) / (0.5 * DELTA_DEMO)) ** 2)
ax[1].fill_between(g * PM / 1000.0, 0, 0.34 * kern, color='#7d3c98', alpha=0.32,
                   lw=0)
ax[1].annotate('a spread of frequencies\nsampled on a steep flank',
               xy=(nu_op * PM / 1000.0, 0.34), xytext=(-0.60, 0.84),
               fontsize=6.5, color='#7d3c98',
               arrowprops=dict(arrowstyle='->', lw=0.7, color='#7d3c98'))
ax[1].axvline(p_true / 1000.0, color='#2980b9', ls=':', lw=0.9)
ax[1].axvline(p_chirp / 1000.0, color='#c0392b', ls=':', lw=0.9)
ax[1].annotate('', xy=(p_chirp / 1000.0, 1.09), xytext=(p_true / 1000.0, 1.09),
               arrowprops=dict(arrowstyle='<|-', lw=1.0, color='k'))
ax[1].text(0.20, 1.14, 'apparent shift %.0f pm' % abs(p_chirp - p_true),
           ha='center', fontsize=6.8)
ax[1].set_ylim(0, 1.28)
ax[1].set_xlim(-0.62, 0.62)
ax[1].set_xlabel('wavelength offset [nm]')
ax[1].set_ylabel('normalised reflectance')
ax[1].set_title('(b) FM-to-AM on the flank', fontsize=9)
ax[1].legend(fontsize=6.8, loc='upper left')

# --- (c) -------------------------------------------------------------------
styles = {0: ('o-', '#c0392b', 'no reference'), 1: ('s-', '#e67e22', '1 reference'),
          2: ('^-', '#2980b9', '2 references'), 3: ('d-', '#28b463', '3 references')}
for n in (0, 1, 2, 3):
    mk, col, lab = styles[n]
    ax[2].plot(ratios, curves[n], mk, color=col, lw=1.2, ms=4, label=lab)
ax[2].axhline(10.0, color='0.3', ls='--', lw=0.8)
ax[2].text(0.97, 0.06, '10 pm target', transform=ax[2].transAxes, fontsize=7,
           color='0.3', ha='right')
ax[2].set_yscale('log')
ax[2].set_xlabel('chirp excursion  ' + r'$\Delta/\mathrm{FWHM}$')
ax[2].set_ylabel('residual Bragg error [pm]')
ax[2].set_title('(c) What the references remove', fontsize=9)
ax[2].legend(fontsize=6.8, loc='upper left')
ax[2].grid(True, which='both', alpha=0.25)

fig.tight_layout()
fig.savefig('figs/fig_s18_source.png', dpi=150, bbox_inches='tight')

print('illustration only: excursion %.2f x FWHM, asymmetry %.2f, apparent shift %.1f pm'
      % (DELTA_DEMO / F, ASYM, p_chirp - p_true))
print('--- residual [pm] vs excursion and reference count ---')
print('  Delta/FWHM ' + ''.join('%9.2f' % r for r in ratios))
for n in (0, 1, 2, 3):
    print('  %d ref%s     ' % (n, ' ' if n == 1 else 's') +
          ''.join('%9.2f' % v for v in curves[n]))
print('saved figs/fig_s18_source.png')
