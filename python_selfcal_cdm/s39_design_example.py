"""s39_design_example.py - one worked instrument, grating by grating.

Table V adds the mechanisms into a single number for a whole array. That is
the right way to compare two designs, and the wrong way to find out which
sensor is the problem. A user who has to certify eight sensors wants to
know what each one is worth, in kelvin, not what their root-sum-square is in
picometres.

So this builds two concrete eight-sensor instruments with every parameter
written down, runs the same serial optical model both times, and reports the
temperature error of each sensor with the mechanisms kept apart.

  naive     what the datasheet suggests: 10 percent gratings at their
            declared 250 pm width, uniform 4 m spacing, wavelengths spread
            evenly over the band, two references inline at the front,
            N = 127 at 25 Mchip/s, no correction.

  designed  every choice taken from the rules of the paper: 1 percent
            gratings narrowed to 100 pm, randomized spacing above 2 m, the
            same wavelength spread, three references on their own stub,
            N = 511 at 100 Mchip/s, sequential deshadowing on.

The optical model is the one used elsewhere: a grating sees the light that
survived two passes through everything in front of it, third-order ghosts land
in whatever bin their delay arithmetic points at, weighted by the triangular
correlation overlap, and the correlation side lobe adds a scaled copy of every
other spectrum. The wavelength axis carries the residual of Section III-D
after the reference fit of s38.

One panel: both arrays side by side per sensor, solid bars for the datasheet
array and hatched bars for the designed one, the net as a cross.

Output: figs/fig_s39_example.pdf, column width.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import common as C
import figstyle as FS

FS.apply(7.0)

VERM, ORAN, GREE, PURP = '#D55E00', '#E69F00', '#009E73', '#CC79A7'
BLUE, GREY = '#0072B2', '0.55'

PM_PER_K = C.TEMP_COEF_PM_PER_C          # 10 pm per kelvin
CA = 4.0 / 3.0 * np.sqrt(2.0 / 3.0)
C_LIGHT, NG = 2.998e8, 1.468
AXIS_PEAK = 6.0                          # smooth axis error before references
# what the reference fit cannot remove, from s15_budget at a chirp excursion
# of 0.2 linewidths: the fit takes out the polynomial part, and this is the
# root-mean-square of what is left over the band
AXIS_FLOOR_PM = {1: 8.57, 2: 4.80, 3: 1.59}

K = 8
NU = np.linspace(-175.0, 175.0, K)       # band positions of the sensors, pm

NAIVE = dict(name='naive', R=0.10, fwhm=250.0, N=127, B=25e6,
             z=4.0 * np.arange(1, K + 1), refs=np.array([-180.0, 180.0]),
             r_ref=0.10, inline=True, peel=False)

DESIGNED = dict(name='designed', R=0.01, fwhm=100.0, N=511, B=100e6,
                z=np.array([2.4, 5.1, 9.3, 12.2, 17.6, 21.0, 27.9, 33.1]),
                refs=np.array([-180.0, 0.0, 180.0]), r_ref=0.01,
                inline=False, peel=True)


def line(nu, nu0, fwhm):
    sig = fwhm / 2.35482
    return np.exp(-0.5 * ((nu - nu0) / sig) ** 2)


def axis_residual(cfg):
    """Wavelength-axis error left at each sensor after the reference fit."""
    sig = cfg['fwhm'] / 2.35482

    def law_a(dl, r):
        return -CA * r * dl * np.exp(-dl ** 2 / (3.0 * sig ** 2))

    def axis(nu):
        return AXIS_PEAK * nu / 200.0

    refs, r_ref = cfg['refs'], cfg['r_ref']
    if cfg['inline']:
        read = np.array([axis(r) + sum(law_a(refs[j] - r, r_ref)
                                       for j in range(i))
                         for i, r in enumerate(refs)])
        seen = np.array([axis(s) + sum(law_a(r - s, r_ref) for r in refs)
                         for s in NU])
    else:
        read = axis(refs)
        seen = axis(NU)
    coef = np.polyfit(refs, read, len(refs) - 1)
    removed = seen - np.polyval(coef, NU)
    # the fit removes the polynomial part of the source error, never all of
    # it. The leading term it cannot touch is one order above the fit, so it
    # is modelled by that Chebyshev polynomial, scaled to the measured floor.
    n = len(refs)
    t = np.cos(n * np.arccos(np.clip(NU / 200.0, -1.0, 1.0)))
    floor = AXIS_FLOOR_PM[n] * t / np.sqrt((t ** 2).mean())
    return removed + floor


def ghost_bins(tau_chips):
    """Third-order arrivals: which sensor bin each path lands in, and which
    three gratings made it."""
    out = {}
    n = len(tau_chips)
    for b in range(n):
        for a in range(n):
            for c in range(n):
                if b < a and b < c:
                    tg = tau_chips[a] - tau_chips[b] + tau_chips[c]
                    for k in range(n):
                        w = 1.0 - abs(tau_chips[k] - tg)
                        if w > 0.02:
                            out.setdefault(k, []).append((a, b, c, w))
    return out


def run(cfg):
    """Per-sensor temperature error, mechanism by mechanism, in kelvin."""
    R, fwhm, N = cfg['R'], cfg['fwhm'], cfg['N']
    nu = np.linspace(-320.0, 320.0, 641)
    tau = 2.0 * NG * cfg['z'] / C_LIGHT * cfg['B']      # delay in chips

    shapes = np.array([line(nu, n0, fwhm) for n0 in NU])
    # serial propagation: two passes through everything in front
    trans = np.ones((K, nu.size))
    for k in range(1, K):
        trans[k] = trans[k - 1] * (1.0 - R * shapes[k - 1]) ** 2
    A = R * shapes * trans

    if cfg['peel']:
        # sequential deshadowing from the source end, as in (13)
        est = np.ones_like(nu)
        A_corr = np.empty_like(A)
        for k in range(K):
            A_corr[k] = A[k] / np.maximum(est, 0.05)
            est = est * (1.0 - np.clip(A_corr[k], 0.0, 0.99)) ** 2
        A_use = A_corr
    else:
        A_use = A

    ghosts = ghost_bins(tau)
    truth = NU.copy()
    err_shadow, err_ghost, err_leak = [], [], []
    for k in range(K):
        base = A_use[k]
        p_sh = C.gauss_fit_peak(nu, base)
        err_shadow.append(p_sh - truth[k])

        g = base.copy()
        for (a, b, c, w) in ghosts.get(k, []):
            prod = shapes[a] * shapes[b] * shapes[c]
            g = g + w * R ** 3 * prod
        p_g = C.gauss_fit_peak(nu, g)
        err_ghost.append(p_g - p_sh)

        leak = -(1.0 / N) * (A_use.sum(axis=0) - A_use[k])
        p_l = C.gauss_fit_peak(nu, g + leak)
        err_leak.append(p_l - p_g)

    err_axis = axis_residual(cfg)
    dz = C_LIGHT / (2.0 * NG * cfg['B'])
    # delay-bin overlap: a rectangular chip gives a triangular correlation
    # peak two chips wide at the base, so a neighbour closer than that leaks
    # a fraction of its own spectrum into this bin before the fit sees it
    err_bin = np.zeros(K)
    for k in range(K):
        add = np.zeros_like(nu)
        for j in range(K):
            if j == k:
                continue
            w = max(0.0, 1.0 - abs(tau[k] - tau[j]))
            if w > 0:
                add = add + w * A_use[j]
        if add.any():
            err_bin[k] = (C.gauss_fit_peak(nu, A_use[k] + add)
                          - C.gauss_fit_peak(nu, A_use[k]))

    parts = dict(shadowing=np.array(err_shadow), ghosts=np.array(err_ghost),
                 leakage=np.array(err_leak), axis=err_axis, delay=err_bin)
    return {k: v / PM_PER_K for k, v in parts.items()}, dz


# ==========================================================================
res_n, dz_n = run(NAIVE)
res_d, dz_d = run(DESIGNED)

ORDER = ['shadowing', 'ghosts', 'leakage', 'axis', 'delay']
COLS = {'shadowing': VERM, 'ghosts': PURP, 'leakage': GREE,
        'axis': BLUE, 'delay': GREY}
WBAR = 0.36
OFF = {'naive': -0.20, 'designed': +0.20}

fig, ax = plt.subplots(figsize=(3.45, 2.55))
x = np.arange(1, K + 1)
for res, cfg in ((res_n, NAIVE), (res_d, DESIGNED)):
    xs = x + OFF[cfg['name']]
    hatched = cfg['name'] == 'designed'
    up = np.zeros(K)
    dn = np.zeros(K)
    for name in ORDER:
        v = res[name]
        pos, neg = np.clip(v, 0, None), np.clip(v, None, 0)
        kw = dict(color=COLS[name], edgecolor='white', lw=0.3)
        if hatched:
            kw = dict(facecolor=COLS[name], alpha=0.55, hatch='/////',
                      edgecolor='white', lw=0.3)
        ax.bar(xs, pos, WBAR, bottom=up, **kw)
        ax.bar(xs, neg, WBAR, bottom=dn, **kw)
        up = up + pos
        dn = dn + neg
    tot = sum(res[n] for n in ORDER)
    ax.plot(xs, tot, 'x', color='0.10', ms=4.6, mew=1.1, zorder=6)

ax.axhline(1.0, color='0.25', ls=(0, (4, 2)), lw=0.85)
ax.axhline(-1.0, color='0.25', ls=(0, (4, 2)), lw=0.85)
ax.axhline(0.0, color='0.6', lw=0.5)
ax.set_xticks(x)
ax.set_xlabel('sensor, in fiber order', labelpad=1.5)
ax.set_ylabel('temperature error [K]', labelpad=1.5)
ax.set_xlim(0.4, K + 0.6)
ax.set_ylim(-2.55, 5.0)
ax.grid(True, axis='y', alpha=0.22)
ax.set_axisbelow(True)
ax.text(0.5, 1.10, '1 K', fontsize=5.6, color='0.25', va='bottom')

handles = [Patch(color=COLS[n], label=n) for n in ORDER]
handles += [Patch(facecolor='0.7', hatch='/////', edgecolor='white',
                  label='hatched: from the rules'),
            Patch(facecolor='0.4', label='solid: from the datasheet'),
            Line2D([], [], ls='none', marker='x', color='0.10', ms=4.6,
                   mew=1.1, label='net error')]
ax.legend(handles=handles, fontsize=5.2, loc='upper left', ncol=2,
          frameon=False, handlelength=1.2, columnspacing=0.8,
          labelspacing=0.2, borderaxespad=0.25, handleheight=0.9)

fig.subplots_adjust(left=0.13, right=0.99, top=0.985, bottom=0.14)
os.makedirs('figs', exist_ok=True)
fig.savefig('figs/fig_s39_example.pdf', bbox_inches='tight', pad_inches=0.01)
fig.savefig('figs/fig_s39_example.png', dpi=200, bbox_inches='tight',
            pad_inches=0.01)
plt.close(fig)

for res, cfg, dz in ((res_n, NAIVE, dz_n), (res_d, DESIGNED, dz_d)):
    tot = sum(res[n] for n in ORDER)
    print('--- %s: R = %.0f%%, FWHM = %.0f pm, N = %d, B = %.0f Mchip/s, '
          'dz = %.2f m ---'
          % (cfg['name'], cfg['R'] * 100, cfg['fwhm'], cfg['N'],
             cfg['B'] / 1e6, dz))
    print('   sensor      ' + ''.join('%7d' % i for i in x))
    for name in ORDER:
        print('   %-11s ' % name + ''.join('%7.2f' % v for v in res[name]))
    print('   %-11s ' % 'net' + ''.join('%7.2f' % v for v in tot))
    print('   worst sensor %.2f K, RMS %.2f K, over 1 K: %d of %d'
          % (np.abs(tot).max(), np.sqrt((tot ** 2).mean()),
             int((np.abs(tot) > 1.0).sum()), K))
    print()
print('saved figs/fig_s39_example.pdf')
