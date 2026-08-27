"""s15_budget.py - the systematic error budget of a low-cost CDM interrogator.

This is the aggregation script: every mechanism studied separately in s3, s7, s8,
s12, s13 and s14 is evaluated for two concrete configurations and put into one
table, so a designer can see which term to attack first.

  BENCH AS PROCURED   3 FBGS gratings, R = 10 percent, 4 m apart, m-sequence
                      N = 127 at 25 Mchip/s, 2 temperature-stabilised references,
                      all gratings inside a +-225 pm window set by the Peltier
                      stages. This is the hardware the project actually owns.

  DESIGNED ARRAY      32 gratings, R = 1 percent, randomised spacing above 2 m,
                      m-sequence N = 511 at 50 Mchip/s, 3 references, sequential
                      deshadowing applied. This is what the design rules of this
                      study recommend once the array grows.

Two rules are kept throughout. First, the source chirp is NOT assumed: the
FM-to-AM term is reported as a function of the chirp excursion expressed in units
of the grating linewidth, because the excursion of the BW10 VCSEL under code
modulation has not been measured yet and belongs on an optical spectrum analyser.
Second, terms are separated into CORRECTABLE systematics (shadowing, wavelength
axis, chirp offset: deterministic, removable by processing or references) and
IRREDUCIBLE ones (noise, code leakage, reference stability), because only the
second group is a genuine floor.
"""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS
FS.apply()

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ                  # 31.25 GHz = 250 pm
SIG = F / 2.35482
C_LIGHT = 2.99792458e8
N_GROUP = 1.468
M_SLOTS = 64
DETUNE_GHZ = 25.0                   # +-200 pm, inside the Peltier tuning range
BAND_HALF = 28.0                    # +-225 pm window of the bench
nu = np.linspace(-2.6 * F, 2.6 * F, M_SLOTS)

# link budget, identical to s12
P_SOURCE_W, LOSS_DB, NEP_W_RTHZ, UNIPOLAR = 1.0e-3, 4.0, 0.5e-12, 0.5


def sigma_noise(nchips, chip_rate):
    p_rx = P_SOURCE_W * 10.0 ** (-LOSS_DB / 10.0)
    return (NEP_W_RTHZ * np.sqrt(chip_rate) / p_rx) / (np.sqrt(nchips) * UNIPOLAR)


# ---------------------------------------------------------------------------
# 1. source chirp, FM-to-AM on the grating edge, after reference correction
# ---------------------------------------------------------------------------
def fmam_residual(chirp_ratio, nref, band_half=BAND_HALF, nsen=8, seed=3):
    """RMS chirp-induced residual [pm] after the references absorb the common part.

    A grating with an asymmetric edge has its reflection peak displaced from the
    nominal Bragg wavelength even with no chirp at all. That displacement is
    static and is removed when the sensor is calibrated, so it is not an
    interrogation error. What this function reports is the CHANGE caused by the
    chirp: each grating is read twice, once with the chirp and once without, and
    the difference is the apparent shift the chirp introduces. The references,
    read the same way, then absorb the smooth, band-dependent part of it, exactly
    as in s3; what survives is the lineshape-dependent remainder.
    """
    rng = np.random.default_rng(seed)
    delta = chirp_ratio * F
    ref_nu = np.linspace(-0.9 * band_half, 0.9 * band_half, max(nref, 1))
    sen_nu = np.sort(rng.uniform(-band_half, band_half, nsen))
    sen_as = rng.uniform(-0.30, 0.30, nsen)
    ref_as = rng.uniform(-0.05, 0.05, max(nref, 1))

    def off(nb):
        # Band dependence of the adiabatic chirp offset. Deliberately NOT a pure
        # polynomial: if it were, a polynomial fit through the references would
        # cancel it exactly and the reference count would look better than it is.
        # The real band dependence is unknown, so a smooth non-polynomial term is
        # included and the fit order is left to work against it.
        x = nb / band_half
        return delta * (0.20 + 0.20 * x + 0.15 * x ** 2 + 0.10 * np.sin(2.5 * x))

    def shift(nb, asym):
        g = np.linspace(nb - 5 * F, nb + 5 * F, 801)
        with_chirp = C.fmam_readout(g, nb, F, delta, mean_off=off(nb), skew=1.2,
                                    shape='tanh', asym=asym)
        without = C.fbg_tanh(g, nb, F, n_side=asym)
        return (C.gauss_fit_peak(g, with_chirp) - C.gauss_fit_peak(g, without)) * PM

    se = np.array([shift(sen_nu[i], sen_as[i]) for i in range(nsen)])
    if nref == 0:
        return float(np.sqrt(np.mean(se ** 2)))
    re = np.array([shift(ref_nu[j], ref_as[j]) for j in range(nref)])
    p = np.polyfit(ref_nu, re, min(nref - 1, 2))
    return float(np.sqrt(np.mean((se - np.polyval(p, sen_nu)) ** 2)))


def lineshape_offset(spread=0.30, nsen=40, seed=13):
    """RMS static displacement [pm] of the reflection peak from the nominal Bragg
    wavelength, caused by edge asymmetry. Reported for context: it is a per-grating
    constant, removed by the sensor calibration, not an interrogation error."""
    rng = np.random.default_rng(seed)
    g = np.linspace(-5 * F, 5 * F, 801)
    e = [C.gauss_fit_peak(g, C.fbg_tanh(g, 0.0, F, n_side=a)) * PM
         for a in rng.uniform(-spread, spread, nsen)]
    return float(np.sqrt(np.mean(np.array(e) ** 2)))


# ---------------------------------------------------------------------------
# 2. nonlinear voltage-to-wavelength axis of the HCG-VCSEL
# ---------------------------------------------------------------------------
SPAN, QUAD, CUBE = 600.0, 0.30, 0.10
_v = np.linspace(0.0, 1.0, 2048)
_nu_true = C.vcsel_sweep(_v, SPAN, quad=QUAD, cube=CUBE)
_nu_lin = _v * SPAN


def axis_residual(nref, band_half=BAND_HALF, centre_frac=0.16, nsen=10, seed=7):
    """RMS residual [pm] of the wavelength axis after an nref-point anchor fit."""
    rng = np.random.default_rng(seed)
    lo, hi = centre_frac * SPAN - band_half, centre_frac * SPAN + band_half
    sen = np.sort(rng.uniform(lo, hi, nsen))
    ref = np.linspace(lo + 0.1 * (hi - lo), hi - 0.1 * (hi - lo), max(nref, 1))
    vst = lambda nb: float(np.interp(nb, _nu_true, _v))
    naive = lambda nb: float(np.interp(vst(nb), _v, _nu_lin))
    raw = np.array([naive(n) for n in sen]) - sen
    if nref == 0:
        return float(np.sqrt(np.mean(raw ** 2)) * PM)
    sv = np.array([vst(n) for n in sen]); rv = np.array([vst(n) for n in ref])
    fit = np.polyval(np.polyfit(rv, ref, min(nref - 1, 2)), sv) - sen
    return float(np.sqrt(np.mean(fit ** 2)) * PM)


# ---------------------------------------------------------------------------
# 3-6. array terms: shadowing, ghosts, code leakage, detector noise
#      (same model as s12, restated here so this script stands alone)
# ---------------------------------------------------------------------------
def array_terms(K, R, nchips, acorr, chip_rate, spacing='random', peel=False, seed=1,
                ntrials=6):
    out = {k: [] for k in ('shadow', 'shadow_corr', 'ghost', 'leak', 'full')}
    for t in range(ntrials):
        rng = np.random.default_rng(seed + t)
        if spacing == 'uniform':
            bins = 1 + max(1, (nchips - 2) // K) * np.arange(K)
        else:
            bins = np.sort(rng.choice(np.arange(1, nchips), size=K, replace=False))
        nub = rng.uniform(-DETUNE_GHZ, DETUNE_GHZ, size=K)
        shapes = np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
        tcum = np.ones((K, M_SLOTS))
        trans = (1.0 - R * shapes) ** 2
        for k in range(1, K):
            tcum[k] = tcum[k - 1] * trans[k - 1]
        clean = R * shapes
        shadowed = R * shapes * tcum
        G = _ghosts(K, nub, R, bins, tcum)
        W = acorr[(bins[:, None] - bins[None, :]) % nchips]; np.fill_diagonal(W, 0.0)
        sg = sigma_noise(nchips, chip_rate)
        fit = lambda S: np.array([(C.gauss_fit_peak(nu, S[k]) - nub[k]) * PM
                                  for k in range(K)])
        rms = lambda e: float(np.sqrt(np.mean(e ** 2)))
        out['shadow'].append(rms(fit(shadowed)))
        out['shadow_corr'].append(rms(_peel(shadowed, nub)))
        out['ghost'].append(rms(fit(clean + G)))
        out['leak'].append(rms(fit(clean + W @ clean + rng.normal(0, sg, clean.shape))))
        S = shadowed + W @ shadowed + G + rng.normal(0, sg, clean.shape)
        out['full'].append(rms(_peel(S, nub) if peel else fit(S)))
    return {k: float(np.mean(v)) for k, v in out.items()}


def _ghosts(K, nub, R, bins, tcum, chunk=40000):
    idx = np.arange(K)
    a, b, c = np.meshgrid(idx, idx, idx, indexing='ij')
    m = (b < a) & (b < c)
    a, b, c = a[m], b[m], c[m]
    gbin = bins[a] - bins[b] + bins[c]
    order = np.argsort(bins)
    pos = np.clip(np.searchsorted(bins[order], gbin), 0, K - 1)
    hit = bins[order][pos] == gbin
    if not hit.any():
        return np.zeros((K, M_SLOTS))
    a, b, c, target = a[hit], b[hit], c[hit], order[pos[hit]]
    cbar = (nub[a] + nub[b] + nub[c]) / 3.0
    spread = (nub[a] - cbar) ** 2 + (nub[b] - cbar) ** 2 + (nub[c] - cbar) ** 2
    amp = R ** 3 * np.exp(-0.5 * spread / SIG ** 2)
    up = np.minimum(np.minimum(a, b), c)
    sig_g = SIG / np.sqrt(3.0)
    out = np.zeros((K, M_SLOTS))
    for i0 in range(0, len(amp), chunk):
        s = slice(i0, i0 + chunk)
        g = (amp[s, None] * tcum[up[s]] *
             np.exp(-0.5 * ((nu[None, :] - cbar[s, None]) / sig_g) ** 2))
        np.add.at(out, target[s], g)
    return out


def _peel(S, nub, floor=0.05):
    K = S.shape[0]
    tcorr = np.ones(M_SLOTS); errs = np.empty(K)
    for k in range(K):
        a, mu, sg, _ = C.gauss_fit_full(nu, S[k] / np.maximum(tcorr, floor))
        errs[k] = (mu - nub[k]) * PM
        tcorr = tcorr * (1.0 - np.clip(a, 0, 0.99) *
                         np.exp(-0.5 * ((nu - mu) / max(sg, 1e-3)) ** 2)) ** 2
    return errs


# ---------------------------------------------------------------------------
# 7. delay-bin resolution
# ---------------------------------------------------------------------------
def resolution_term(K, chip_rate, d_min, nchips, ntrials=10, seed=5):
    side = -1.0 / nchips
    errs = []
    for t in range(ntrials):
        rng = np.random.default_rng(seed + t)
        z = np.cumsum(np.concatenate([[2.0], d_min * rng.uniform(1.0, 2.5, K - 1)]))
        tc = (2.0 * N_GROUP * z / C_LIGHT) * chip_rate
        nub = rng.uniform(-DETUNE_GHZ, DETUNE_GHZ, size=K)
        A = np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
        for k in range(K):
            dt = np.round(tc[k]) - tc
            w = np.where(np.abs(dt) < 1.0, np.clip(1 - np.abs(dt), 0, None),
                         side)
            errs.append((C.gauss_fit_peak(nu, w @ A) - nub[k]) * PM)
    return float(np.sqrt(np.mean(np.array(errs) ** 2)))


# ---------------------------------------------------------------------------
# configurations
# ---------------------------------------------------------------------------
_m7 = 1.0 - 2.0 * C._mls01(7); AC7 = C.periodic_xcorr(_m7, _m7)
_m9 = 1.0 - 2.0 * C._mls01(9); AC9 = C.periodic_xcorr(_m9, _m9)

CHIRPS = [0.1, 0.2, 0.3]
LINESHAPE_OFF = None   # filled after the helpers are defined
REF_STAB_PM = C.TEMP_STAB_C * C.TEMP_COEF_PM_PER_C     # 0.1 C x 10 pm/C

CFG = {
    'bench as procured': dict(K=3, R=0.10, nchips=127, acorr=AC7, chip_rate=25e6,
                              d_min=4.0, nref=2, peel=False, spacing='uniform'),
    'designed array': dict(K=32, R=0.01, nchips=511, acorr=AC9, chip_rate=100e6,
                           d_min=2.0, nref=3, peel=True, spacing='random'),
}

ROWS = ['source chirp, FM-AM residual', 'wavelength axis nonlinearity',
        'spectral shadowing', 'multiple reflections (ghosts)',
        'code leakage', 'detector noise', 'delay-bin resolution',
        'reference temperature stability']
CORRECTABLE = {'source chirp, FM-AM residual', 'wavelength axis nonlinearity',
               'spectral shadowing'}

LINESHAPE_OFF = lineshape_offset()

budget = {}
for name, cfg in CFG.items():
    at = array_terms(cfg['K'], cfg['R'], cfg['nchips'], cfg['acorr'],
                     cfg['chip_rate'], spacing=cfg['spacing'], peel=cfg['peel'])
    sg = sigma_noise(cfg['nchips'], cfg['chip_rate'])
    rng = np.random.default_rng(99)
    ne = [(C.gauss_fit_peak(nu, cfg['R'] * np.exp(-0.5 * ((nu - x) / SIG) ** 2) +
                            rng.normal(0, sg, M_SLOTS)) - x) * PM
          for x in rng.uniform(-DETUNE_GHZ, DETUNE_GHZ, 200)]
    budget[name] = {
        'source chirp, FM-AM residual': {c: fmam_residual(c, cfg['nref']) for c in CHIRPS},
        'wavelength axis nonlinearity': axis_residual(cfg['nref']),
        'spectral shadowing': at['shadow_corr'] if cfg['peel'] else at['shadow'],
        'multiple reflections (ghosts)': at['ghost'],
        'code leakage': at['leak'],
        'detector noise': float(np.std(ne)),
        'delay-bin resolution': resolution_term(cfg['K'], cfg['chip_rate'],
                                                cfg['d_min'], cfg['nchips']),
        'reference temperature stability': REF_STAB_PM,
        '_full': at['full'], '_peel': cfg['peel'],
    }

# ---------------------------------------------------------------------------
# printed budget
# ---------------------------------------------------------------------------
print('SYSTEMATIC ERROR BUDGET, low-cost CDM interrogation of an FBG array')
print('gratings: FWHM %.0f pm, detuning spread +-%.0f pm; estimator: Gaussian fit'
      % (C.FBG_FWHM_PM, DETUNE_GHZ * PM))
for name, cfg in CFG.items():
    print()
    print('=== %s ===' % name.upper())
    print('    K=%d, R=%.0f%%, code m-sequence N=%d, %.0f Mchip/s, min spacing %.1f m,'
          % (cfg['K'], cfg['R'] * 100, cfg['nchips'], cfg['chip_rate'] / 1e6, cfg['d_min']))
    print('    %d references, %s spacing, deshadowing %s'
          % (cfg['nref'], cfg['spacing'], 'on' if cfg['peel'] else 'off'))
    b = budget[name]
    print('    %-34s %10s   %s' % ('term', 'RMS [pm]', 'character'))
    tot_sq = 0.0
    for r in ROWS:
        v = b[r]
        if isinstance(v, dict):
            txt = ' / '.join('%.2f' % v[c] for c in CHIRPS)
            print('    %-34s %10s   correctable, shown for chirp/FWHM = %s'
                  % (r, txt, ' / '.join(str(c) for c in CHIRPS)))
            tot_sq += v[CHIRPS[1]] ** 2
        else:
            tag = 'correctable systematic' if r in CORRECTABLE else 'irreducible'
            txt = '<0.01' if v < 0.01 else '%.3f' % v
            print('    %-34s %10s   %s' % (r, txt, tag))
            tot_sq += v ** 2
    print('    %-34s %10.2f   root-sum-square, chirp/FWHM = %.1f'
          % ('TOTAL', np.sqrt(tot_sq), CHIRPS[1]))
    irr = np.sqrt(sum(b[r] ** 2 for r in ROWS if r not in CORRECTABLE
                      and not isinstance(b[r], dict)))
    print('    %-34s %10.2f   what stays after every correction'
          % ('irreducible part', irr))
    print('    %-34s %10.2f   array terms only, end to end'
          % ('cross-check, array model', b['_full']))
    print('    %-34s %10.2f   static, removed by sensor calibration'
          % ('(grating lineshape offset)', LINESHAPE_OFF))

# ---------------------------------------------------------------------------
# how many temperature-stabilised references are worth having
# ---------------------------------------------------------------------------
print()
print('=== REFERENCE COUNT: chirp residual [pm] over the +-225 pm bench window ===')
print('   references   fit order   chirp/FWHM 0.1   0.2   0.3   0.5')
for nr, order in [(0, None), (1, 'constant'), (2, 'linear'), (3, 'quadratic'),
                  (4, 'quadratic')]:
    vals = [fmam_residual(c, nr) for c in (0.1, 0.2, 0.3, 0.5)]
    print('   %10d   %9s   %14.2f  %5.2f  %5.2f  %5.2f'
          % (nr, order if order else 'none', *vals))
print('   the references remove the smooth band dependence of the chirp offset;')
print('   what stays scales with the chirp and with the lineshape spread')

# ---------------------------------------------------------------------------
# figure: horizontal bars per term, both configurations
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.1, 2.9))
y = np.arange(len(ROWS))
h = 0.38
vals = {}
for name in CFG:
    vals[name] = np.array([budget[name][r][CHIRPS[1]] if isinstance(budget[name][r], dict)
                           else budget[name][r] for r in ROWS])
FLOOR = 0.012
ax.barh(y + h / 2, np.maximum(vals['bench as procured'], FLOOR), height=h,
        color='#c0392b', label='bench as procured (K=3, R=10%, N=127, 2 refs)')
ax.barh(y - h / 2, np.maximum(vals['designed array'], FLOOR), height=h,
        color='#2980b9', label='designed array (K=32, R=1%, N=511, 3 refs, deshadowed)')
for name, dy, col in [('bench as procured', h / 2, '#c0392b'),
                      ('designed array', -h / 2, '#2980b9')]:
    for yi, v in zip(y, vals[name]):
        if v < 0.01:
            ax.text(FLOOR * 1.4, yi + dy, 'below 0.01', va='center', fontsize=6.5, color=col)
ax.set_yticks(y); ax.set_yticklabels(ROWS, fontsize=6.4)
ax.set_xscale('log'); ax.set_xlim(0.01, 100); ax.set_ylim(-0.8, len(ROWS) - 0.2)
ax.axvline(10.0, color='0.3', ls='--', lw=0.9)
ax.text(10.6, len(ROWS) - 0.6, '10 pm target', fontsize=7.5, color='0.3')
ax.set_xlabel('RMS contribution to the Bragg-wavelength error [pm]')
ax.set_title('Systematic error budget, chirp excursion = %.1f x linewidth' % CHIRPS[1], fontsize=7)
ax.grid(True, axis='x', which='both', alpha=0.25)
ax.legend(fontsize=5.8, loc='lower right', frameon=False)
fig.subplots_adjust(left=0.26, right=0.98, top=0.91, bottom=0.13); plt.savefig('figs/fig_s15_budget.png', dpi=140)
print()
print('saved figs/fig_s15_budget.png')
