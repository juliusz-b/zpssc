"""s23_theory.py - closed-form first-order laws for the coded readout, verified.

The simulations so far produced scaling laws numerically. This script derives
them, and then checks the derivations against the model to the percent level.
The common starting point is a linear response argument: every estimator used
here maps a spectrum to a number, and for a small additive perturbation p(nu) on
a Gaussian line A*g(nu), g = exp(-nu^2/2 s^2), the reading shifts by a linear
functional of p. For the two standard estimators the functional is explicit:

  centroid:  delta = (1/(A s sqrt(2 pi))) * Int nu p(nu) dnu
  LSQ fit :  delta = (2/(A s sqrt(pi)))  * Int (nu/s) p(nu) g(nu) dnu
             (unweighted least squares, free amplitude/centre/width/baseline;
              on a symmetric window the centre couples only to the odd part)

LAW A, pairwise shadowing. A downstream grating read through one upstream
grating of reflectance R detuned by D has p = -2 R A g(nu) g(nu - D), and the
integrals close:

  centroid:  delta(D) = -(R D / sqrt(2)) exp(-D^2 / 4 s^2)
  LSQ     :  delta(D) = -(4/3) sqrt(2/3) R D exp(-D^2 / 3 s^2)

Both vanish for co-tuned and for far-detuned neighbours; the worst case sits at
D* = s sqrt(2) (centroid) or s sqrt(3/2) (LSQ), with peak magnitude 0.607 R s
and 0.809 R s. The dangerous neighbour is the one detuned by about one line
sigma; a co-tuned neighbour biases nothing.

LAW B, mean multiple-access bias. An m-sequence correlator leaks every other
grating into the wanted bin with weight -1/N. Averaging the leaked background
over detunings uniform on (-W, W) and applying the LSQ functional gives

  delta_bar(nu_k) = ((K-1) s^2 / (W N)) *
                    [ exp(-(W-nu_k)^2/4s^2) - exp(-(W+nu_k)^2/4s^2) ]

an odd, monotone pull toward the outside of the band, linear in (K-1)/N and
independent of reflectivity. For the centroid the same average is exactly
(K-1) nu_k / N: a pure stretch of the wavelength axis by 1 + (K-1)/N.

COROLLARY. Because the mean bias is an odd, nearly linear function of position,
the two reference gratings that calibrate the wavelength axis remove most of it
as a by-product. Verified below: 9.4 pm of leakage bias falls to 2.3 pm at
K = 32, N = 127.

INEQUALITY C, capacity-refresh product. Leakage caps the per-band array at
K <= eps N / c_L, and ADC-limited equivalent-time sampling sets the refresh rate
f_r = f_s / (n_s M N). The code length cancels from the product:

  K * f_r <= eps f_s / (c_L n_s M)

so gratings-per-band times refresh rate is an invariant of the receiver, set by
the converter rate and the accuracy target, not by N or the chip rate.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS
FS.apply()


PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
W = 25.0
NCH = 127
_MS = 1.0 - 2.0 * C._mls01(7)
ACORR = C.periodic_xcorr(_MS, _MS)

CL = (4.0 / 3.0) * np.sqrt(2.0 / 3.0)          # LSQ prefactor, ~1.089


def centroid_full(x, y):
    return float(np.sum(x * y) / np.sum(y))


def lsq_full(x, y, mu0):
    f = lambda t, a, m, s, c: a * np.exp(-0.5 * ((t - m) / s) ** 2) + c
    try:
        p, _ = curve_fit(f, x, y, p0=[max(y.max(), 1e-9), mu0, SIG, 0.0],
                         maxfev=20000)
        return float(p[1])
    except Exception:
        return mu0


# ---------------------------------------------------------------------------
# LAW A verification
# ---------------------------------------------------------------------------
g = np.linspace(-5 * F, 5 * F, 3000)
R_A = 0.01
Dfr = np.linspace(0.1, 3.0, 25)
simc, siml = [], []
for Df in Dfr:
    D = Df * SIG
    meas = np.exp(-0.5 * (g / SIG) ** 2) * \
        (1 - R_A * np.exp(-0.5 * ((g - D) / SIG) ** 2)) ** 2
    simc.append(centroid_full(g, meas) * PM)
    siml.append(lsq_full(g, meas, 0.0) * PM)
simc, siml = np.array(simc), np.array(siml)
thc = -R_A * Dfr * SIG / np.sqrt(2.0) * np.exp(-(Dfr * SIG) ** 2 / (4 * SIG ** 2)) * PM
thl = -CL * R_A * Dfr * SIG * np.exp(-(Dfr * SIG) ** 2 / (3 * SIG ** 2)) * PM

err_c = float(np.max(np.abs(simc / thc - 1)))
err_l = float(np.max(np.abs(siml / thl - 1)))

# ---------------------------------------------------------------------------
# LAW B verification and the reference corollary
# ---------------------------------------------------------------------------
K = 32
R_B = 0.10
M = 256
nu = np.linspace(-2.6 * F, 2.6 * F, M)


def leak_trial(rng, refs=False):
    bins = np.sort(rng.choice(np.arange(1, NCH), size=K, replace=False))
    nub = rng.uniform(-W, W, size=K)
    if refs:
        nub[0], nub[1] = -0.9 * W, 0.9 * W
    A = R_B * np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
    Wl = ACORR[(bins[:, None] - bins[None, :]) % NCH]
    np.fill_diagonal(Wl, 0.0)
    S = A + Wl @ A
    mus = np.array([lsq_full(nu, S[k], nub[k]) for k in range(K)])
    return nub, mus


pos, err = [], []
for s_ in range(30):
    nub, mus = leak_trial(np.random.default_rng(s_))
    pos.extend(nub); err.extend((mus - nub) * PM)
pos, err = np.array(pos), np.array(err)
rms_raw = float(np.sqrt(np.mean(err ** 2)))

res_ref = []
for s_ in range(30):
    nub, mus = leak_trial(np.random.default_rng(1000 + s_), refs=True)
    p = np.polyfit(mus[:2], nub[:2], 1)
    cor = np.polyval(p, mus[2:]) - nub[2:]
    res_ref.append(np.sqrt(np.mean((cor * PM) ** 2)))
rms_ref = float(np.mean(res_ref))

edges = np.linspace(-W, W, 11)
cent = 0.5 * (edges[:-1] + edges[1:])
binned = np.array([err[(pos >= a) & (pos < b)].mean()
                   for a, b in zip(edges[:-1], edges[1:])])
pref = (K - 1) * SIG ** 2 / (W * NCH)
lawB = pref * (np.exp(-(W - cent) ** 2 / (4 * SIG ** 2)) -
               np.exp(-(W + cent) ** 2 / (4 * SIG ** 2))) * PM
nu_fine = np.linspace(-W, W, 300)
lawB_fine = pref * (np.exp(-(W - nu_fine) ** 2 / (4 * SIG ** 2)) -
                    np.exp(-(W + nu_fine) ** 2 / (4 * SIG ** 2))) * PM

# ---------------------------------------------------------------------------
# INEQUALITY C check against the design table
# ---------------------------------------------------------------------------
EPS, FS, NSAMP, MST, C_L = 10.0, 20e6, 4, 64, 27.0
bound = EPS * FS / (C_L * NSAMP * MST)
table = {63: (24.6, 1033.0), 127: (49.2, 513.0), 255: (94.2, 306.0),
         511: (191.6, 127.0)}
products = {N: k * f for N, (k, f) in table.items()}

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(7.1, 2.4))

ax[0].plot(Dfr, thc, '-', color='#2980b9', lw=1.3, label='centroid, closed form')
ax[0].plot(Dfr[::2], simc[::2], 'o', color='#2980b9', ms=4, mfc='none',
           label='centroid, simulated')
ax[0].plot(Dfr, thl, '-', color='#c0392b', lw=1.3, label='LSQ fit, closed form')
ax[0].plot(Dfr[::2], siml[::2], 's', color='#c0392b', ms=4, mfc='none',
           label='LSQ fit, simulated')
ax[0].axvline(np.sqrt(1.5), color='#c0392b', ls=':', lw=0.8)
ax[0].axvline(np.sqrt(2.0), color='#2980b9', ls=':', lw=0.8)
ax[0].text(1.36, -0.14, r'$\Delta^{*}$', fontsize=8, color='0.3')
ax[0].set_xlabel(r'neighbour detuning  $\Delta/\sigma$')
ax[0].set_ylabel('pairwise shadowing bias [pm]  (R = 1%)')
ax[0].set_title('(a) Law A, exact to first order', fontsize=9)
ax[0].legend(fontsize=6.6, loc='lower right'); ax[0].grid(True, alpha=0.25)

ax[1].plot(pos, err, '.', ms=2, color='0.75', alpha=0.5)
ax[1].plot(cent, binned, 'o', color='#2980b9', ms=5, label='binned mean, simulated')
ax[1].plot(nu_fine, lawB_fine, '-', color='#c0392b', lw=1.4, label='Law B, closed form')
ax[1].set_xlabel(r'grating position in the band  $\nu_k$ [GHz]')
ax[1].set_ylabel('multiple-access bias [pm]')
ax[1].set_title('(b) Law B: an odd pull, K = 32, N = 127', fontsize=9)
ax[1].text(0.03, 0.84, 'after 2-reference axis fit:\n%.1f pm -> %.1f pm'
           % (rms_raw, rms_ref), transform=ax[1].transAxes, fontsize=7,
           color='#28b463')
ax[1].legend(fontsize=6.8, loc='lower right'); ax[1].grid(True, alpha=0.25)

Ns = sorted(products)
ax[2].plot(Ns, [products[N] / 1e3 for N in Ns], 'o-', color='#2980b9',
           label=r'simulated  $K_{\max} \cdot f_r$')
ax[2].axhline(bound / 1e3, color='#c0392b', ls='--', lw=1.2,
              label=r'bound  $\epsilon f_s / (c_L n_s M)$')
ax[2].set_xscale('log')
ax[2].set_xticks(Ns); ax[2].set_xticklabels([str(N) for N in Ns], fontsize=7.5)
ax[2].minorticks_off()
ax[2].set_ylim(0, bound / 1e3 * 1.25)
ax[2].set_xlabel('code length N')
ax[2].set_ylabel(r'capacity-refresh product [10$^3$ sensor$\cdot$Hz]')
ax[2].set_title('(c) Inequality C: N cancels', fontsize=9)
ax[2].legend(fontsize=6.8, loc='lower right'); ax[2].grid(True, which='both', alpha=0.25)

fig.tight_layout()
fig.savefig('figs/fig_s23_theory.png', dpi=150, bbox_inches='tight')

# ---------------------------------------------------------------------------
# printed verification
# ---------------------------------------------------------------------------
print('=== LAW A: pairwise shadowing bias, R = %.0f%% ===' % (R_A * 100))
print('  centroid: worst deviation from closed form %.1f%%' % (100 * err_c))
print('  LSQ fit : worst deviation from closed form %.1f%%' % (100 * err_l))
print('  worst-case bias: 0.607 R sigma (centroid), 0.809 R sigma (LSQ)')
print('  at R = 10%%, FWHM = 250 pm: %.1f pm at a neighbour detuning of %.0f pm'
      % (0.809 * 0.10 * SIG * PM, np.sqrt(1.5) * SIG * PM))
print()
print('=== LAW B: mean multiple-access bias, K = %d, N = %d ===' % (K, NCH))
print('   bin centre [GHz]   simulated [pm]   closed form [pm]')
for c_, b_, t_ in zip(cent, binned, lawB):
    print('   %14.1f   %14.2f   %16.2f' % (c_, b_, t_))
print('  RMS raw %.2f pm; after the 2-reference axis fit %.2f pm (%.1fx)'
      % (rms_raw, rms_ref, rms_raw / rms_ref))
print()
print('=== INEQUALITY C: capacity-refresh product ===')
print('  bound: eps*fs/(c_L*n_s*M) = %.0f sensor*Hz per band' % bound)
for N in Ns:
    print('   N = %4d: K_max * f_r = %.0f  (%.0f%% of the bound)'
          % (N, products[N], 100 * products[N] / bound))
print('saved figs/fig_s23_theory.png')
