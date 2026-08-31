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

One panel, one axis pair. The eight datasheet curves are drawn thin with a
shaded band spanning their spread, the designed curves lie on the diagonal.

Output: figs/fig_s40_temperature.pdf, column width.
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

PM_PER_K = C.TEMP_COEF_PM_PER_C
K = 8
NU0 = np.linspace(-175.0, 175.0, K)      # nominal band positions, pm
SWEEP = np.linspace(-200.0, 200.0, 41)   # excursion of the moving sensor, pm

CFG = {
    'datasheet': dict(R=0.10, fwhm=250.0, N=127, peel=False),
    'designed': dict(R=0.01, fwhm=100.0, N=511, peel=True),
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


tk = SWEEP / PM_PER_K
rec = {}
for name, cfg in CFG.items():
    rows = []
    for k in range(K):
        rec_pm = []
        for d in SWEEP:
            centres = NU0.copy()
            centres[k] = NU0[k] + d
            rec_pm.append(recover(centres, cfg, k))
        rows.append((np.array(rec_pm) - NU0[k]) / PM_PER_K)
    rec[name] = np.array(rows)

fig, ax = plt.subplots(figsize=(3.45, 3.30))

arr = rec['datasheet']
ax.fill_between(tk, arr.min(axis=0), arr.max(axis=0), color='#D55E00',
                alpha=0.16, lw=0, zorder=1)
cmap = plt.get_cmap('Oranges')
for k in range(K):
    ax.plot(tk, arr[k], lw=0.75, color=cmap(0.45 + 0.5 * k / (K - 1.0)),
            zorder=2)
for k in range(K):
    ax.plot(tk, rec['designed'][k], lw=0.9, color='#0072B2', zorder=3)
ax.plot([-20, 20], [-20, 20], color='0.25', lw=0.8, ls=(0, (4, 2)), zorder=4)

# the worst point, named in kelvin
dep = arr - tk[None, :]
ks, i = np.unravel_index(np.abs(dep).argmax(), dep.shape)
ax.annotate('up to %.1f K\noff the diagonal' % abs(dep[ks, i]),
            xy=(tk[i], arr[ks, i]), xytext=(1.5, -13.0), fontsize=5.6,
            color='#D55E00', ha='left', va='center',
            arrowprops=dict(arrowstyle='-|>', color='#D55E00', lw=0.7,
                            mutation_scale=7))

ax.set_xlim(-21, 21)
ax.set_ylim(-21, 21)
ax.set_aspect('equal')
ax.set_xlabel('true change $\\Delta T$ [K]', labelpad=1.5)
ax.set_ylabel('recovered change $\\Delta\\hat T$ [K]', labelpad=1.5)
ax.grid(True, alpha=0.22)
ax.set_axisbelow(True)

handles = [
    Line2D([], [], color=cmap(0.7), lw=1.2,
           label='datasheet, one line per sensor'),
    Patch(facecolor='#D55E00', alpha=0.16, label='spread of the eight'),
    Line2D([], [], color='#0072B2', lw=1.2,
           label='designed, all eight'),
    Line2D([], [], color='0.25', lw=0.8, ls=(0, (4, 2)),
           label='recovered $=$ true'),
]
ax.legend(handles=handles, fontsize=5.3, loc='upper left', frameon=False,
          handlelength=1.5, labelspacing=0.22, borderaxespad=0.25)

fig.subplots_adjust(left=0.135, right=0.99, top=0.99, bottom=0.115)
os.makedirs('figs', exist_ok=True)
fig.savefig('figs/fig_s40_temperature.pdf', bbox_inches='tight',
            pad_inches=0.01)
fig.savefig('figs/fig_s40_temperature.png', dpi=200, bbox_inches='tight',
            pad_inches=0.01)
plt.close(fig)

for name in CFG:
    d = rec[name] - tk[None, :]
    j = np.abs(d).max(axis=0).argmax()
    print('%-10s worst departure %5.2f K (at a true change of %+.1f K), '
          'RMS %5.2f K, span across sensors at the worst point %5.2f K'
          % (name, np.abs(d).max(), tk[j], np.sqrt((d ** 2).mean()),
             (d.max(axis=0) - d.min(axis=0)).max()))
print('saved figs/fig_s40_temperature.pdf')
