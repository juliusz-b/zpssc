"""s40_temperature.py - recovered temperature against the real one.

A calibration curve is the plot a sensor user actually reads. It should be a
straight line, and the interesting question is where it stops being one.

Each sensor is swept across the band in turn while the other seven stay at
their nominal wavelengths, which is the case of one hot spot on a fiber that
is otherwise in equilibrium. As the moving grating passes a neighbour, their
detuning sweeps through the range where Law A is strongest, so the recovered
temperature leaves the line and comes back. The wobble is not an offset. An
offset would calibrate away, and this does not, because its size depends on
where the sensor happens to be.

One column-wide figure. The upper panel is the calibration curve itself for
both arrays of the worked example, one line per sensor. The lower panel,
sharing the temperature axis, is the departure from the diagonal, where the
kelvins are legible.

Output: figs/fig_s40_temperature.pdf, column width.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import common as C
import figstyle as FS

FS.apply(7.0)

PM_PER_K = C.TEMP_COEF_PM_PER_C
K = 8
NU0 = np.linspace(-175.0, 175.0, K)      # nominal band positions, pm
SWEEP = np.linspace(-200.0, 200.0, 41)   # excursion of the moving sensor, pm

CFG = {
    'datasheet': dict(R=0.10, fwhm=250.0, N=127, peel=False, cmap='Oranges'),
    'designed': dict(R=0.01, fwhm=100.0, N=511, peel=True, cmap='Blues'),
}


def line(nu, nu0, fwhm):
    sig = fwhm / 2.35482
    return np.exp(-0.5 * ((nu - nu0) / sig) ** 2)


def recover(centres, cfg, which):
    """Fitted Bragg position of one grating, given where all of them sit."""
    R, fwhm, N = cfg['R'], cfg['fwhm'], cfg['N']
    nu = np.linspace(-560.0, 560.0, 1121)
    shapes = np.array([line(nu, c, fwhm) for c in centres])
    trans = np.ones((K, nu.size))
    for k in range(1, K):
        trans[k] = trans[k - 1] * (1.0 - R * shapes[k - 1]) ** 2
    A = R * shapes * trans

    if cfg['peel']:
        est = np.ones_like(nu)
        corr = np.empty_like(A)
        for k in range(K):
            corr[k] = A[k] / np.maximum(est, 0.05)
            est = est * (1.0 - np.clip(corr[k], 0.0, 0.99)) ** 2
        A = corr

    leak = -(1.0 / N) * (A.sum(axis=0) - A[which])
    return C.gauss_fit_peak(nu, A[which] + leak)


fig, (top, bot) = plt.subplots(2, 1, figsize=(3.45, 3.55), sharex=True,
                               gridspec_kw=dict(height_ratios=[2.0, 1.0],
                                                hspace=0.08))
store = {}
tk = SWEEP / PM_PER_K

for name, cfg in CFG.items():
    cmap = plt.get_cmap(cfg['cmap'])
    for k in range(K):
        rec_pm = []
        for d in SWEEP:
            centres = NU0.copy()
            centres[k] = NU0[k] + d
            rec_pm.append(recover(centres, cfg, k))
        rec_k = (np.array(rec_pm) - NU0[k]) / PM_PER_K
        col = cmap(0.45 + 0.5 * k / (K - 1.0))
        top.plot(tk, rec_k, lw=0.9, color=col)
        bot.plot(tk, rec_k - tk, lw=0.9, color=col)
        store.setdefault(name, []).append(rec_k - tk)

top.plot([-20, 20], [-20, 20], color='0.25', lw=0.8, ls=(0, (4, 2)), zorder=0)
top.set_xlim(-21, 21)
top.set_ylim(-21, 21)
top.set_ylabel('recovered change [K]', labelpad=1.5)
top.grid(True, alpha=0.22)
top.set_axisbelow(True)
handles = [Line2D([], [], color=plt.get_cmap('Oranges')(0.7), lw=1.2,
                  label='datasheet array, one line per sensor'),
           Line2D([], [], color=plt.get_cmap('Blues')(0.7), lw=1.2,
                  label='designed array, eight lines on the diagonal'),
           Line2D([], [], color='0.25', lw=0.8, ls=(0, (4, 2)),
                  label='recovered $=$ true')]
top.legend(handles=handles, fontsize=5.4, loc='upper left', frameon=False,
           handlelength=1.5, labelspacing=0.2, borderaxespad=0.3)

bot.axhline(0.0, color='0.25', lw=0.8, ls=(0, (4, 2)))
for lim in (1.0, -1.0):
    bot.axhline(lim, color='0.55', lw=0.7, ls=(0, (1.5, 2)))
bot.text(-20.3, 1.15, '1 K', fontsize=5.4, color='0.45', va='bottom')
bot.set_xlabel('true change [K]', labelpad=1.5)
bot.set_ylabel('recovered $-$ true [K]', labelpad=1.5)
bot.set_ylim(-3.6, 5.2)
bot.set_yticks([-2, 0, 2, 4])
bot.grid(True, alpha=0.22)
bot.set_axisbelow(True)

fig.subplots_adjust(left=0.135, right=0.99, top=0.99, bottom=0.10)
os.makedirs('figs', exist_ok=True)
fig.savefig('figs/fig_s40_temperature.pdf', bbox_inches='tight',
            pad_inches=0.01)
fig.savefig('figs/fig_s40_temperature.png', dpi=200, bbox_inches='tight',
            pad_inches=0.01)
plt.close(fig)

for name in CFG:
    arr = np.array(store[name])
    i = np.abs(arr).max(axis=0).argmax()
    print('%-10s worst departure %5.2f K (at a true change of %+.1f K), '
          'RMS %5.2f K, span across sensors at the worst point %5.2f K'
          % (name, np.abs(arr).max(), tk[i], np.sqrt((arr ** 2).mean()),
             (arr.max(axis=0) - arr.min(axis=0)).max()))
print('saved figs/fig_s40_temperature.pdf')
