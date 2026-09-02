"""s41_spectra.py - the reconstructed line, before and after deshadowing.

Everything else in the paper reports a fitted number. This shows the object
the number is fitted to, at three depths in the array.

The first grating is unshadowed and its measurement is already the true line,
which is what lets the recursion of (13) start. By the eighth the line has
lost a third of its height and leans to one side, and the correction puts
back both.

The last panel answers a question the first three raise. Spacing decides
whether third-order arrivals collide with a grating, so it ought to change
the reconstructed line, and it does not. A ghost carries R cubed against a
line carrying R. Where they pile up, at the far end of the array, their
sum reaches 7.1 percent of the line on the uniform grid and 1.8 percent
on the randomized one, a factor of four that leaves the drawn line
looking the same.
Ghost placement is worth doing, but the payoff shows up in the fitted error
over many layouts, in Figs. 4 and 6, not in any single spectrum.

Output: figs/fig_s41_spectra.pdf, full text width.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import common as C
import figstyle as FS

FS.apply(7.0)

VERM, GREE, BLUE = '#D55E00', '#009E73', '#0072B2'

K = 8
R = 0.10
FWHM = 250.0
N = 127
B = 100e6                      # 1.02 m resolution, so the delays are resolved
C_LIGHT, NG = 2.998e8, 1.468
NU0 = np.linspace(-175.0, 175.0, K)
SHOW = [0, 3, 7]

LAYOUT = {
    'uniform, 4 m': 4.0 * np.arange(1, K + 1),
    'randomized': np.array([2.4, 5.1, 9.3, 12.2, 17.6, 21.0, 27.9, 33.1]),
}


def line(nu, nu0):
    sig = FWHM / 2.35482
    return np.exp(-0.5 * ((nu - nu0) / sig) ** 2)


def build(z):
    """Measured and corrected spectra, and the ghost term on its own."""
    nu = np.linspace(-520.0, 520.0, 1041)
    tau = 2.0 * NG * z / C_LIGHT * B
    shapes = np.array([line(nu, c) for c in NU0])

    trans = np.ones((K, nu.size))
    for k in range(1, K):
        trans[k] = trans[k - 1] * (1.0 - R * shapes[k - 1]) ** 2
    A = R * shapes * trans

    ghost = np.zeros_like(A)
    for b in range(K):
        for a in range(K):
            for c in range(K):
                if b < a and b < c:
                    tg = tau[a] - tau[b] + tau[c]
                    for k in range(K):
                        w = 1.0 - abs(tau[k] - tg)
                        if w > 0.02:
                            ghost[k] += (w * R ** 3 * shapes[a] * shapes[b]
                                         * shapes[c])

    meas = np.empty_like(A)
    for k in range(K):
        meas[k] = A[k] + ghost[k] - (1.0 / N) * (A.sum(axis=0) - A[k])

    est = np.ones_like(nu)
    corr = np.empty_like(A)
    for k in range(K):
        corr[k] = meas[k] / np.maximum(est, 0.05)
        est = est * (1.0 - np.clip(corr[k], 0.0, 0.99)) ** 2
    return nu, meas, corr, ghost


nu, meas, corr, ghost_u = build(LAYOUT['uniform, 4 m'])
_, _, _, ghost_r = build(LAYOUT['randomized'])
lam = nu / 1000.0

fig, ax = plt.subplots(1, 4, figsize=(7.16, 1.95))
for col, k in enumerate(SHOW):
    a = ax[col]
    true = R * line(nu, NU0[k])
    a.plot(lam, true, color='0.55', lw=2.6, alpha=0.55, label='true $R_k(\\lambda)$')
    a.plot(lam, meas[k], color=VERM, lw=1.2, label='measured $S_k$')
    a.plot(lam, corr[k], color=GREE, lw=1.1, ls=(0, (3.5, 2)),
           label='deshadowed $\\widehat S_k$')
    a.axvline(NU0[k] / 1000.0, color='0.3', lw=0.6, ls=(0, (1.5, 2.5)))
    loss = 100.0 * (1.0 - meas[k].max() / true.max())
    a.text(-0.5, 0.1255, 'grating %d' % (k + 1), fontsize=6.2, color='0.2',
           va='top')
    a.text(-0.5, 0.1135, '%.0f%% of the height gone' % loss, fontsize=5.6,
           color='0.4', va='top')
    a.set_xlim(-0.52, 0.52)
    a.set_ylim(0, 0.132)
    a.set_yticks([0, 0.05, 0.10])
    a.set_xlabel('wavelength offset [nm]', labelpad=1.5)
    a.grid(False, alpha=0.2)
    a.set_axisbelow(True)
    if col:
        a.set_yticklabels([])
ax[0].set_ylabel('reflectance', labelpad=1.5)
ax[0].legend(fontsize=5.5, loc='center right', frameon=True,
             handlelength=1.5, labelspacing=0.18, borderaxespad=0.22)

# --- the ghost term on its own, at the scale it actually has --------------
a = ax[3]
true8 = R * line(nu, NU0[7])
a.plot(lam, true8, color='0.55', lw=2.6, alpha=0.55,
       label='the line $A_8(\\lambda)$')
a.plot(lam, ghost_u[7], color=VERM, lw=1.2,
       label='$+$ summed $P_g$, uniform')
a.plot(lam, ghost_r[7], color=BLUE, lw=1.2,
       label='$+$ summed $P_g$, randomized')
a.set_yscale('log')
a.set_ylim(1e-8, 0.6)
a.set_xlim(-0.52, 0.52)
a.set_xlabel('wavelength offset [nm]', labelpad=1.5)
a.set_ylabel('reflectance', labelpad=1.5)
a.text(-0.5, 0.30, 'grating 8: what ghosts add on top', fontsize=6.2,
       color='0.2', va='top',
       bbox=dict(fc='white', ec='none', pad=0.6, alpha=0.8))
a.text(-0.5, 0.105, 'drawn alone, log scale',
       fontsize=5.2,
       color='0.4', va='top',
       bbox=dict(fc='white', ec='none', pad=0.6, alpha=0.8))
leg = a.legend(fontsize=5.2, loc='lower center', frameon=True,
               framealpha=0.85, edgecolor='none', handlelength=1.3,
               labelspacing=0.18, borderaxespad=0.2)
leg.set_zorder(7)
a.grid(False, which='both', alpha=0.2)
a.set_axisbelow(True)

fig.subplots_adjust(left=0.062, right=0.995, top=0.975, bottom=0.205,
                    wspace=0.16)
os.makedirs('figs', exist_ok=True)
fig.savefig('figs/fig_s41_spectra.pdf', bbox_inches='tight', pad_inches=0.01)
fig.savefig('figs/fig_s41_spectra.png', dpi=200, bbox_inches='tight',
            pad_inches=0.01)
plt.close(fig)

for k in SHOW:
    true = R * line(nu, NU0[k])
    m = true > 0.02 * true.max()
    print('grating %d: height left %.0f%%, measured %+6.2f pm, '
          'corrected %+6.2f pm, shape mismatch after correction %.2f%%'
          % (k + 1, 100 * meas[k].max() / true.max(),
             C.gauss_fit_peak(nu, meas[k]) - NU0[k],
             C.gauss_fit_peak(nu, corr[k]) - NU0[k],
             100 * np.sqrt(np.mean((corr[k][m] - true[m]) ** 2)) / true.max()))
print('ghost peak at grating 8: uniform %.2e, randomized %.2e, ratio %.1f'
      % (ghost_u[7].max(), ghost_r[7].max(),
         ghost_u[7].max() / ghost_r[7].max()))
print('that peak is %.2f percent of the line'
      % (100 * ghost_u[7].max() / (R * line(nu, NU0[7])).max()))
print('saved figs/fig_s41_spectra.pdf')
