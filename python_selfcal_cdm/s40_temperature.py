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

The two panels use the two configurations of the worked example. At the
datasheet reflectivity the curve departs by kelvins and every sensor departs
differently, since each has a different set of neighbours in front of it. At
one percent the same eight curves lie on the line to within the width of the
stroke.

Output: figs/fig_s40_temperature.pdf, full text width.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import common as C
import figstyle as FS

FS.apply(7.0)

PM_PER_K = C.TEMP_COEF_PM_PER_C
K = 8
NU0 = np.linspace(-175.0, 175.0, K)      # nominal band positions, pm
SWEEP = np.linspace(-200.0, 200.0, 41)   # excursion of the moving sensor, pm

CFG = {
    'datasheet': dict(R=0.10, fwhm=250.0, N=127, peel=False, col='#D55E00'),
    'designed': dict(R=0.01, fwhm=100.0, N=511, peel=True, col='#0072B2'),
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


fig, ax = plt.subplots(1, 3, figsize=(7.16, 2.28),
                       gridspec_kw=dict(width_ratios=[1.0, 1.0, 1.12]))
cmap = plt.get_cmap('viridis')
store = {}

for a, (name, cfg) in zip(ax[:2], CFG.items()):
    for k in range(K):
        true_pm, rec_pm = [], []
        for d in SWEEP:
            centres = NU0.copy()
            centres[k] = NU0[k] + d
            true_pm.append(centres[k])
            rec_pm.append(recover(centres, cfg, k))
        true_k = (np.array(true_pm) - NU0[k]) / PM_PER_K
        rec_k = (np.array(rec_pm) - NU0[k]) / PM_PER_K
        a.plot(true_k, rec_k, lw=1.0, color=cmap(k / (K - 1.0)))
        store.setdefault(name, []).append(rec_k - true_k)
    a.plot([-20, 20], [-20, 20], color='0.25', lw=0.8, ls=(0, (4, 2)),
           zorder=0)
    a.set_xlim(-21, 21)
    a.set_ylim(-21, 21)
    a.set_aspect('equal')
    a.set_xlabel('true change [K]', labelpad=1.5)
    a.set_title('(%s) %s' % ('a' if name == 'datasheet' else 'b', name),
                fontsize=7.4)
    a.grid(True, alpha=0.22)
    a.set_axisbelow(True)

ax[0].set_ylabel('recovered change [K]', labelpad=1.5)
ax[0].text(-19.5, 19.0,
           'one line per sensor,' '\n' 'dashed is recovered $=$ true',
           fontsize=5.6, color='0.3', va='top')

# --- (c) the departure, where it is legible -------------------------------
tk = (SWEEP) / PM_PER_K
for name, cfg in CFG.items():
    arr = np.array(store[name])
    ax[2].fill_between(tk, arr.min(axis=0), arr.max(axis=0), alpha=0.25,
                       color=cfg['col'], lw=0)
    ax[2].plot(tk, arr.mean(axis=0), color=cfg['col'], lw=1.3, label=name)
ax[2].axhline(0.0, color='0.25', lw=0.8, ls=(0, (4, 2)))
for lim in (1.0, -1.0):
    ax[2].axhline(lim, color='0.55', lw=0.7, ls=(0, (1.5, 2)))
ax[2].text(-19.5, 1.15, '1 K', fontsize=5.6, color='0.45', va='bottom')
ax[2].set_xlim(-21, 21)
ax[2].set_xlabel('true change [K]', labelpad=1.5)
ax[2].set_ylabel('recovered $-$ true [K]', labelpad=1.5)
ax[2].set_title('(c) the departure, all eight', fontsize=7.4)
ax[2].legend(fontsize=5.6, loc='lower left', frameon=False, handlelength=1.5,
             labelspacing=0.2, borderaxespad=0.3)
ax[2].grid(True, alpha=0.22)
ax[2].set_axisbelow(True)

fig.subplots_adjust(left=0.062, right=0.995, top=0.885, bottom=0.175,
                    wspace=0.30)
os.makedirs('figs', exist_ok=True)
fig.savefig('figs/fig_s40_temperature.pdf', bbox_inches='tight',
            pad_inches=0.01)
fig.savefig('figs/fig_s40_temperature.png', dpi=200, bbox_inches='tight',
            pad_inches=0.01)
plt.close(fig)

for name in CFG:
    arr = np.array(store[name])
    print('%-10s worst departure %5.2f K, RMS %5.2f K, span across sensors '
          'at the worst point %5.2f K'
          % (name, np.abs(arr).max(), np.sqrt((arr ** 2).mean()),
             (arr.max(axis=0) - arr.min(axis=0)).max()))
print('saved figs/fig_s40_temperature.pdf')
