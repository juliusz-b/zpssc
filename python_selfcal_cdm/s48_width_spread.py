"""s48_width_spread.py - what happens to Law A and to the array budget when
the gratings do not share one linewidth.

Fabricated gratings never have exactly equal FWHM, so this script answers the
reviewer question the equal-width Law A leaves open. Three results.

  (a) The width-aware generalization of Law A. Repeating the first-order
      projection with unequal widths sigma_j (neighbor) and sigma_k (read
      grating) gives

        bias = -4 sqrt(2) R dl sigma_j sigma_k^2 / (2 sigma_j^2 + sigma_k^2)^(3/2)
               * exp(-dl^2 / (2 sigma_j^2 + sigma_k^2))

      which collapses to the equal-width law for sigma_j = sigma_k (the
      prefactor becomes (4/3) sqrt(2/3) and the exponent 3 sigma^2). Checked
      against an unrestricted least-squares fit at R = 0.5%: agreement to
      0.3%. Panel (a) draws both at R = 10% for three width ratios.

  (b) Worst-case bias versus width ratio, three conventions: the first-order
      law, the full two-pass model with an unrestricted fit (the convention
      of s26 and of the 9.0-vs-8.6 pm sentence in the text), and the windowed
      peak fit the acquisition protocol actually uses (gauss_fit_full). The
      window amplifies the worst case by 15-45%, most for narrow neighbors,
      because a narrow notch bites the core of the line and the window fits
      exactly that core. The most dangerous neighbor is narrower than the
      read grating by about sqrt(2), not a wider one.

  (c) The eight-grating array with exact serial shadowing and m-sequence code
      leakage, FWHM of every grating drawn uniformly within a fabrication
      spread of 0-30%. Raw error at R = 10%, the same after sequential
      deshadowing, and raw error at R = 1%. The budget moves by <10%, and
      the deshadowing correction is unaffected, because it divides the
      spectra point by point and assumes neither a width nor a shape.

Output: figs/fig_s48_widths.pdf/.png, drawn at 7.1 in for a full-width figure.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import common as C
import figstyle as FS
FS.apply()

SIG = 250.0 / 2.35482                    # pm, read grating
NU = np.linspace(-5 * SIG, 5 * SIG, 1800)


def gauss(x, mu, sig):
    return np.exp(-0.5 * ((x - mu) / sig) ** 2)


def law(R, dl, sj, sk=SIG):
    s2 = 2 * sj ** 2 + sk ** 2
    return -4 * np.sqrt(2) * R * dl * sj * sk ** 2 / s2 ** 1.5 * np.exp(-dl ** 2 / s2)


def _g4(x, A, mu, s, b):
    return A * np.exp(-0.5 * ((x - mu) / s) ** 2) + b


def fit_unrestricted(y):
    p, _ = curve_fit(_g4, NU, y, p0=[y.max() - y.min(), NU[np.argmax(y)], SIG,
                                     y.min()], maxfev=20000)
    return p[1]


def model(R, d, sj):
    return gauss(NU, 0, SIG) * (1 - R * gauss(NU, d, sj)) ** 2


R = 0.10
DLS = np.linspace(10, 460, 46)

fig, ax = plt.subplots(1, 3, figsize=(7.1, 2.35))

# --- (a) bias curves: full model (solid) vs law (dashed) --------------------
cols = {0.7: FS.VERM, 1.0: FS.BLUE, 1.3: FS.GREEN}
for ratio in (0.7, 1.0, 1.3):
    sj = SIG * ratio
    num = np.array([fit_unrestricted(model(R, d, sj)) for d in DLS])
    ax[0].plot(DLS, num, color=cols[ratio], lw=1.4,
               label=r'$\sigma_j/\sigma_k=%.1f$' % ratio)
    ax[0].plot(DLS, law(R, DLS, sj), color=cols[ratio], lw=1.0, ls='--')
ax[0].axhline(0, color='0.6', lw=0.6)
ax[0].set_xlabel(r'$\Delta\lambda_{jk}$ (pm)')
ax[0].set_ylabel(r'$\delta\lambda_{k\leftarrow j}$ (pm)')
ax[0].legend(fontsize=6.2, loc='lower right', handlelength=1.6)

# --- (b) worst case vs width ratio, three conventions -----------------------
ratios = np.linspace(0.6, 1.4, 17)
mx_law, mx_unr, mx_win = [], [], []
for r_ in ratios:
    sj = SIG * r_
    mx_law.append(np.max(np.abs(law(R, DLS, sj))))
    mx_unr.append(np.max(np.abs([fit_unrestricted(model(R, d, sj)) for d in DLS])))
    mx_win.append(np.max(np.abs([C.gauss_fit_full(NU, model(R, d, sj))[1]
                                 for d in DLS])))
ax[1].plot(ratios, mx_win, 'o-', color=FS.PURPLE, ms=2.8, label='windowed peak fit')
ax[1].plot(ratios, mx_unr, 's-', color=FS.BLUE, ms=2.8, label='full model, unrestricted fit')
ax[1].plot(ratios, mx_law, '--', color=FS.VERM, label='first-order law')
ax[1].set_xlabel(r'width ratio $\sigma_j/\sigma_k$')
ax[1].set_ylabel(r'worst-case $|\delta\lambda_{k\leftarrow j}|$ (pm)')
ax[1].legend(fontsize=6.2, handlelength=1.6)

# --- (c) K=8 array, FWHM spread ---------------------------------------------
NCH = 127
MS = 1.0 - 2.0 * C._mls01(7)
ACORR = C.periodic_xcorr(MS, MS)
MSTEP = 64
GRID = np.linspace(-512.0, 512.0, MSTEP)


def run_array(K, Rr, spread, rng, deshadow):
    fw = 250.0 * (1 + rng.uniform(-spread, spread, K))
    sg = fw / 2.35482
    nub = rng.uniform(-200, 200, K)
    bins = np.sort(rng.choice(np.arange(1, NCH), size=K, replace=False))
    shapes = np.array([gauss(GRID, nub[k], sg[k]) for k in range(K)])
    tcum = np.ones((K, MSTEP))
    for k in range(1, K):
        tcum[k] = tcum[k - 1] * (1.0 - Rr * shapes[k - 1]) ** 2
    A = Rr * shapes * tcum
    W = ACORR[(bins[:, None] - bins[None, :]) % NCH]
    np.fill_diagonal(W, 0.0)
    S = A + W @ A
    if deshadow:
        T = np.ones(MSTEP)
        Sc = np.zeros_like(S)
        for k in range(K):
            Sc[k] = S[k] / T
            T = T * np.clip(1.0 - Sc[k], 0.05, None) ** 2
        S = Sc
    errs = [C.gauss_fit_full(GRID, S[k])[1] - nub[k] for k in range(K)]
    return float(np.sqrt(np.mean(np.square(errs))))


NT = 30
spreads = [0.0, 0.10, 0.20, 0.30]
for label, Rr, ds, col, mk in (('$R=10\\%$, raw', 0.10, False, FS.VERM, 'o'),
                               ('$R=10\\%$, deshadowed', 0.10, True, FS.BLUE, 's'),
                               ('$R=1\\%$, raw', 0.01, False, FS.GREEN, '^')):
    row = [np.mean([run_array(8, Rr, sp, np.random.default_rng(900 + t), ds)
                    for t in range(NT)]) for sp in spreads]
    ax[2].plot([100 * s for s in spreads], row, mk + '-', color=col, ms=3.4,
               label=label)
    print(label, ['%.2f' % v for v in row])
ax[2].axhline(10, color='0.4', ls='--', lw=0.8)
ax[2].text(29, 10.6, '10 pm', fontsize=6.2, color='0.4', ha='right')
ax[2].set_yscale('log')
ax[2].set_xlabel('FWHM spread (%)')
ax[2].set_ylabel('array RMS error (pm)')
ax[2].legend(fontsize=6.2, handlelength=1.6)

for a, letter in zip(ax, 'abc'):
    a.text(0.02, 1.06, letter, transform=a.transAxes, fontsize=9,
           fontweight='bold', va='bottom')

fig.tight_layout(pad=0.4, w_pad=1.2)
fig.savefig('figs/fig_s48_widths.png', dpi=300, bbox_inches='tight', pad_inches=0.01)
fig.savefig('figs/fig_s48_widths.pdf', bbox_inches='tight', pad_inches=0.01)
print('saved figs/fig_s48_widths.pdf')
