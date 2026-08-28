"""s24_concept.py - end-to-end measurement chain used as Fig. 1.

The figure follows one quantity through the interrogator. At every nominal
wavelength step the VCSEL transmits one period of the same binary sequence.
The photodiode records a sum of delayed, wavelength-weighted copies. A
periodic correlation converts that row into peaks at the grating delays.
Repeating the operation over the sweep makes a delay-wavelength map. A cut at
one delay is one grating spectrum, whose fitted centre is the measurement.

The visual encoding is deliberately restrained. Blue denotes the shared
measurement chain. Orange denotes the selected second grating and follows it
from the fiber through the delay peak, map cut, and fitted Bragg wavelength.
No colour is used as a wavelength scale.
"""
import os
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np

import common as C
import figstyle as FS

warnings.filterwarnings('ignore')
FS.apply(base=7.1)

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
NCH = 127
M = 40
K = 4
R = 0.05
BINS = np.array([17, 46, 77, 107])
NU_B = np.array([-16.0, 9.0, -5.0, 19.0])
SELECTED = 1

MSEQ = 1.0 - 2.0 * C._mls01(7)
MSEQ = np.roll(MSEQ, -int(np.argmax(MSEQ == 1.0)))
ACORR = C.periodic_xcorr(MSEQ, MSEQ)


def build_measurement():
    """Return the simulated record, correlated rows, and selected fit."""
    nu = np.linspace(-1.9 * F, 1.9 * F, M)
    lam_pm = nu * PM
    shapes = np.exp(-0.5 * ((nu[None, :] - NU_B[:, None]) / SIG) ** 2)

    transmission = np.ones((K, M))
    for k in range(1, K):
        transmission[k] = transmission[k - 1] * \
            (1.0 - R * shapes[k - 1]) ** 2
    primary = R * shapes * transmission

    weights = ACORR[(np.arange(NCH)[:, None] - BINS[None, :]) % NCH]
    correlated = (weights @ primary).T
    correlated /= R

    osr = 2
    code01 = np.repeat(0.5 * (1.0 + MSEQ), osr)
    components = np.empty((M, K, NCH * osr))
    for m in range(M):
        for k in range(K):
            components[m, k] = primary[k, m] * np.roll(code01, BINS[k] * osr)
    rng = np.random.default_rng(7)
    record = components.sum(axis=1)
    record += rng.normal(0.0, 0.0025, record.shape)

    m_sel = int(np.argmin(np.abs(nu - NU_B[SELECTED])))
    spectrum = correlated[:, BINS[SELECTED]]
    amp, mu, sig, base = C.gauss_fit_full(nu, spectrum)
    nu_fit = np.linspace(nu.min(), nu.max(), 400)
    fit = base + amp * np.exp(-0.5 * ((nu_fit - mu) / sig) ** 2)
    return nu, lam_pm, primary, record, components, correlated, m_sel, \
        spectrum, nu_fit, fit, mu


nu, lam_pm, primary, record, components, corr, m_sel, spec, nu_fit, fit, mu = \
    build_measurement()


def title(ax, letter, text):
    ax.set_title(r'$\bf{%s}$  %s' % (letter, text), loc='left', pad=4,
                 fontsize=7.3)


def arrow(ax, xy1, xy2, color='0.35', lw=0.9, mutation=8, connection='arc3'):
    ax.add_patch(FancyArrowPatch(xy1, xy2, arrowstyle='-|>', color=color,
                                 lw=lw, mutation_scale=mutation,
                                 connectionstyle=connection))


fig = plt.figure(figsize=(7.1, 3.12))
outer = fig.add_gridspec(2, 1, height_ratios=[0.86, 1.18], hspace=0.46)
top = outer[0].subgridspec(1, 2, width_ratios=[1.15, 1.85], wspace=0.20)
bottom = outer[1].subgridspec(1, 4, width_ratios=[1.15, 1.00, 1.43, 0.72],
                              wspace=0.38)

# ---------------------------------------------------------------------------
# (a) One code period at each wavelength step
# ---------------------------------------------------------------------------
ax = fig.add_subplot(top[0])
title(ax, 'a', 'One code period per wavelength')
ax.set_xlim(-0.18, 4.30)
ax.set_ylim(-0.35, 2.95)
ax.axis('off')

bits = (0.5 * (1.0 + MSEQ[:10])).astype(int)
step_y = [0.0, 0.72, 1.44, 2.16]
step_labels = [r'$\lambda_1$', r'$\lambda_2$', r'$\lambda_3$', r'$\lambda_M$']
for s, (y0, lab) in enumerate(zip(step_y, step_labels)):
    x0 = 0.05 + 1.02 * s
    for q, bit in enumerate(bits):
        ax.add_patch(Rectangle((x0 + q * 0.092, y0), 0.088, 0.30,
                               facecolor=FS.BLUE if bit else 'white',
                               edgecolor=FS.BLUE, linewidth=0.42))
    ax.text(x0 + 0.45, y0 + 0.40, lab, ha='center', va='bottom',
            fontsize=7.0, color=FS.BLUE)
    if s < 3:
        ax.plot([x0 + 0.92, x0 + 1.00], [y0 + 0.15, step_y[s + 1] + 0.15],
                color='0.55', lw=0.55, ls=(0, (2, 1.5)))
ax.annotate('', xy=(-0.10, 2.70), xytext=(-0.10, -0.08),
            arrowprops=dict(arrowstyle='-|>', color='0.35', lw=0.7))
ax.text(-0.16, 1.32, 'nominal\nwavelength', ha='right', va='center',
        fontsize=6.4, color='0.35')
ax.annotate('', xy=(4.20, -0.20), xytext=(0.02, -0.20),
            arrowprops=dict(arrowstyle='-|>', color='0.35', lw=0.7))
ax.text(4.20, -0.28, 'time', ha='right', va='top', fontsize=6.4, color='0.35')
ax.text(2.02, 2.76, r'filled/open cells = $c(t)=1/0$', ha='center',
        fontsize=6.4, color='0.30')

# ---------------------------------------------------------------------------
# (b) Fiber and delayed returns
# ---------------------------------------------------------------------------
ax = fig.add_subplot(top[1])
title(ax, 'b', 'Return delay identifies the FBG')
ax.set_xlim(0, 10)
ax.set_ylim(-0.3, 4.6)
ax.axis('off')

source = FancyBboxPatch((0.05, 2.29), 1.45, 0.94,
                        boxstyle='round,pad=0.06,rounding_size=0.08',
                        fc='#EAF3F8', ec=FS.BLUE, lw=0.9)
ax.add_patch(source)
ax.text(0.775, 2.76, 'swept\nVCSEL', ha='center', va='center', fontsize=6.8)
ax.add_patch(Circle((2.15, 2.76), 0.25, fc='white', ec='0.35', lw=0.8))
ax.text(2.15, 2.76, 'C', ha='center', va='center', fontsize=6.2, color='0.35')
arrow(ax, (1.52, 2.76), (1.88, 2.76), FS.BLUE)
ax.text(1.70, 2.98, r'$\lambda_m,\ c(t)$', ha='center', va='bottom',
        fontsize=6.2, color=FS.BLUE)
ax.plot([2.40, 9.75], [2.76, 2.76], color='0.30', lw=1.45)

z = [3.35, 4.85, 6.45, 8.55]
cols = ['0.48', FS.ORANGE, '0.48', '0.48']
for k, (zk, col) in enumerate(zip(z, cols), start=1):
    for dx in (-0.055, 0.0, 0.055):
        ax.plot([zk + dx, zk + dx], [2.50, 3.02], color=col, lw=1.0)
    ax.text(zk, 3.16, 'FBG %d' % k, ha='center', fontsize=6.2,
            color=col)
    ax.text(zk, 2.30, r'$z_%d$' % k, ha='center', fontsize=6.4,
            color=col)

ax.plot([2.15, 2.15], [2.51, 1.10], color='0.35', lw=0.9)
pd = FancyBboxPatch((1.55, 0.42), 1.20, 0.67,
                    boxstyle='round,pad=0.04,rounding_size=0.06',
                    fc='white', ec='0.35', lw=0.8)
ax.add_patch(pd)
ax.text(2.15, 0.755, 'photodiode', ha='center', va='center', fontsize=6.7)
ax.text(5.95, 4.38, r'at $\lambda_m$: amplitude $A_k(\lambda_m)$',
        ha='center', fontsize=6.5, color='0.30')
for k, zk in enumerate(z):
    col = FS.ORANGE if k == SELECTED else FS.BLUE
    arrow(ax, (zk - 0.10, 3.95 - 0.17 * k),
          (2.43, 3.95 - 0.17 * k), col, lw=0.72, mutation=6)
ax.text(6.15, 1.03,
        r'$p_m(t)=\sum_k A_k(\lambda_m)c(t-\tau_k)+n(t)$',
        ha='center', va='center', fontsize=7.2, color='0.18')
ax.text(6.15, 0.42, r'$\tau_k=2n_gz_k/c$  $\longrightarrow$  position',
        ha='center', fontsize=6.8, color=FS.ORANGE)

# ---------------------------------------------------------------------------
# (c) Raw detector row
# ---------------------------------------------------------------------------
ax_raw = fig.add_subplot(bottom[0])
title(ax_raw, 'c', r'Raw PD sum')
t = np.arange(record.shape[1]) / 2.0
scale = record[m_sel].max()
for k in range(K):
    ax_raw.plot(t, components[m_sel, k] / scale, color='0.72', lw=0.42,
                alpha=0.75)
ax_raw.plot(t, record[m_sel] / scale, color='0.15', lw=0.75,
            label=r'recorded $p_m(t)$')
ax_raw.set_xlim(0, NCH)
ax_raw.set_ylim(-0.10, 1.20)
ax_raw.set_xticks([0, 64, 127])
ax_raw.set_yticks([])
ax_raw.set_xlabel('time [chips]')
ax_raw.text(0.03, 0.90, 'sum seen by the PD', transform=ax_raw.transAxes,
            fontsize=6.3, color='0.20')
ax_raw.text(0.03, 0.06, 'grey: individual echoes', transform=ax_raw.transAxes,
            fontsize=5.9, color='0.45')
ax_raw.spines['left'].set_visible(False)

# ---------------------------------------------------------------------------
# (d) Correlation at one wavelength
# ---------------------------------------------------------------------------
ax_cor = fig.add_subplot(bottom[1])
title(ax_cor, 'd', r'Correlation $\rightarrow$ delay peaks')
delay = np.arange(NCH)
row = corr[m_sel]
ax_cor.plot(delay, row, color=FS.BLUE, lw=0.85)
for k, bk in enumerate(BINS):
    col = FS.ORANGE if k == SELECTED else FS.BLUE
    ax_cor.vlines(bk, 0, row[bk], color=col, lw=1.0)
    ax_cor.plot(bk, row[bk], 'o', ms=3.2, color=col, mec='white', mew=0.35)
    ax_cor.text(bk, row[bk] + 0.055, r'$\tau_%d$' % (k + 1), ha='center',
                fontsize=6.3, color=col)
ax_cor.set_xlim(0, NCH)
ax_cor.set_ylim(min(-0.025, row.min() * 1.25), row.max() * 1.28)
ax_cor.set_xticks([0, 64, 127])
ax_cor.set_yticks([])
ax_cor.set_xlabel(r'delay $\tau$ [chips]')
ax_cor.text(0.04, 0.08, r'peak height $=A_k(\lambda_m)$',
            transform=ax_cor.transAxes, fontsize=6.0, color='0.30')
ax_cor.spines['left'].set_visible(False)

# ---------------------------------------------------------------------------
# (e) Delay-wavelength map
# ---------------------------------------------------------------------------
ax_map = fig.add_subplot(bottom[2])
title(ax_map, 'e', r'Rows $\rightarrow$ map')
display = np.clip(corr, 0.0, None)
ax_map.imshow(display, origin='lower', aspect='auto', cmap='Blues',
              extent=[0, NCH, lam_pm[0], lam_pm[-1]],
              vmin=0.0, vmax=np.percentile(display, 99.8))
for k, bk in enumerate(BINS):
    col = FS.ORANGE if k == SELECTED else '0.22'
    ax_map.text(bk, lam_pm[-1] - 22, 'FBG %d' % (k + 1), ha='center',
                va='top', fontsize=5.9, color=col)
ax_map.axvline(BINS[SELECTED], color=FS.ORANGE, lw=1.0, ls=(0, (3, 2)))
ax_map.set_xlim(5, 118)
ax_map.set_xlabel(r'delay $\tau$  $\longrightarrow$  grating position')
ax_map.set_ylabel(r'$\lambda_m-\lambda_0$ [pm]')
ax_map.tick_params(pad=1)
ax_map.text(0.04, 0.06, r'cut at $\tau_2$', transform=ax_map.transAxes,
            fontsize=6.4, color=FS.ORANGE,
            bbox=dict(boxstyle='round,pad=0.18', fc='white', ec=FS.ORANGE,
                      lw=0.6, alpha=0.92))

# ---------------------------------------------------------------------------
# (f) One column is one spectrum
# ---------------------------------------------------------------------------
ax_fit = fig.add_subplot(bottom[3], sharey=ax_map)
title(ax_fit, 'f', r'Cut $\rightarrow$ spectrum')
den = max(spec.max(), 1e-12)
ax_fit.plot(spec / den, lam_pm, 'o', color=FS.ORANGE, ms=2.6,
            mec='white', mew=0.3)
ax_fit.plot(fit / den, nu_fit * PM, color=FS.ORANGE, lw=1.0)
lam_b_pm = mu * PM
ax_fit.axhline(lam_b_pm, color='0.25', lw=0.75, ls=(0, (3, 2)))
ax_fit.text(0.04, lam_b_pm + 25, r'$\hat{\lambda}_{B,2}$', fontsize=7.3,
            color='0.18', va='bottom')
ax_fit.set_xlim(-0.05, 1.08)
ax_fit.set_xticks([0, 1])
ax_fit.set_xlabel('reflectivity')
ax_fit.tick_params(axis='y', labelleft=False)
ax_fit.spines['left'].set_visible(False)

os.makedirs('figs', exist_ok=True)
fig.savefig('figs/fig_s24_concept.pdf', bbox_inches='tight', pad_inches=0.025)
fig.savefig('figs/fig_s24_concept.png', dpi=300, bbox_inches='tight',
            pad_inches=0.025)
plt.close(fig)

print('selected FBG 2: true %.1f pm, fitted %.1f pm' %
      (NU_B[SELECTED] * PM, lam_b_pm))
print('saved figs/fig_s24_concept.pdf and figs/fig_s24_concept.png')
