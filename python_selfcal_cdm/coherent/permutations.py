"""All 40320 orders of the eight gratings of Fig. 3(b) (R = 10 %, FWHM 175..325 pm, co-tuned), on deterministic
coupled-mode direct spectra (s70_field), without noise.

For every order: the largest Omega_k of (eq:omega) over the array, whether every return has a single maximum
at zero detuning, and the largest shift of a return maximum when all gratings upstream of the read grating move
together by +5 pm (the quantity bounded by (eq:globalgain) times 5 pm).

python permutations.py -> figs/fig_order_permutations.pdf/.png, numbers on stdout,
                          cache/order_perm.npz (read by fig_validation.py for Fig. 3c)
"""
import sys, itertools, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
import s70_field as F
import figstyle as FS

FIGS = ROOT.parent / 'figs'
OUT = ROOT / 'cache' / 'order_perm.npz'
R0, DELTA = 0.10, 5.0
widths = np.linspace(175.0, 325.0, 8)
xf = np.linspace(-350.0, 350.0, 1401)
dx = xf[1] - xf[0]
length, kappa, _, _ = F.calibrate('bl250', R0)


def refl(w, shift):
    return abs(F.grating_r(F.LAM0 + (xf - shift) * 1e-12, F.LAM0, length * 250 / w, kappa * w / 250, sections=160)) ** 2


def maxima(A):
    """Number of local maxima above 30 % of the peak and the refined position of the largest one, per row."""
    up = A[:, 1:-1] >= A[:, :-2]
    dn = A[:, 1:-1] >= A[:, 2:]
    hi = A[:, 1:-1] > 0.3 * A.max(axis=1, keepdims=True)
    n = (up & dn & hi).sum(axis=1)
    i = np.argmax(A, axis=1)
    i = np.clip(i, 1, A.shape[1] - 2)
    r = np.arange(A.shape[0])
    ym, y0, yp = A[r, i - 1], A[r, i], A[r, i + 1]
    den = ym - 2 * y0 + yp
    off = np.where(den != 0, 0.5 * (ym - yp) / np.where(den != 0, den, 1), 0.0)
    return n, xf[i] + off * dx


def main():
    t0 = time.time()
    q0 = np.array([refl(w, 0.0) for w in widths])          # (8, nx) own spectra, no shift
    q5 = np.array([refl(w, DELTA) for w in widths])        # upstream spectra shifted by +5 pm
    T0 = (1 - q0) ** 2
    T5 = (1 - q5) ** 2
    sig = widths / 2.354820045
    gamma = 1 / sig ** 2
    p = 2 * R0 / (1 - R0) * gamma
    perms = np.array(list(itertools.permutations(range(8))), dtype=np.int8)   # (40320, 8)
    P = len(perms)
    omega_max = np.empty(P)
    n_split = np.zeros(P, dtype=int)        # gratings with more than one maximum at delta = 0
    shift_max = np.empty(P)                 # largest |offset of the maximum from lambda_B| over the array for delta = +5 pm
    shift_last = np.empty(P)
    for b in range(0, P, 1024):
        pp = perms[b:b + 1024]
        n = len(pp)
        cp = np.cumsum(p[pp], axis=1)
        om = np.concatenate([np.zeros((n, 1)), cp[:, :-1]], axis=1) / gamma[pp]
        omega_max[b:b + n] = om.max(axis=1)
        prod0 = np.cumprod(T0[pp], axis=1)                # (n, 8, nx)
        prod5 = np.cumprod(T5[pp], axis=1)
        one = np.ones((n, 1, len(xf)))
        A0 = q0[pp] * np.concatenate([one, prod0[:, :-1]], axis=1)
        A5 = q0[pp] * np.concatenate([one, prod5[:, :-1]], axis=1)
        nm0, x0 = maxima(A0.reshape(-1, len(xf)))
        nm5, x5 = maxima(A5.reshape(-1, len(xf)))
        nm0 = nm0.reshape(n, 8); x0 = x0.reshape(n, 8); x5 = x5.reshape(n, 8)
        n_split[b:b + n] = (nm0 > 1).sum(axis=1)
        sh = np.abs(x5)                                   # offset of the return maximum from lambda_B
        shift_max[b:b + n] = sh.max(axis=1)
        shift_last[b:b + n] = sh[:, -1]
    print('computed %d orders in %.0f s' % (P, time.time() - t0))
    np.savez(OUT, perms=perms, omega_max=omega_max, n_split=n_split, shift_max=shift_max, shift_last=shift_last,
             widths=widths, R0=R0, delta=DELTA)

    i_wide = int(np.where((perms == np.arange(8)[::-1]).all(axis=1))[0][0])
    i_narrow = int(np.where((perms == np.arange(8)).all(axis=1))[0][0])
    i_best = int(np.argmin(omega_max))
    ok = omega_max < 1
    print('orders with Omega_k < 1 for every grating: %d of %d (%.1f %%)' % (ok.sum(), P, 100 * ok.mean()))
    print('orders with every return single: %d (%.1f %%)' % ((n_split == 0).sum(), 100 * (n_split == 0).mean()))
    print('single returns among Omega<1: %d of %d, split among Omega>=1: %d of %d' % (
        (n_split[ok] == 0).sum(), ok.sum(), (n_split[~ok] > 0).sum(), (~ok).sum()))
    print('widest first: Omega_max %.2f, shift %.1f pm, bound %.1f pm' % (
        omega_max[i_wide], shift_max[i_wide], DELTA * omega_max[i_wide] / (1 - omega_max[i_wide])))
    print('narrowest first: Omega_max %.2f, shift %.1f pm, split gratings %d' % (omega_max[i_narrow], shift_max[i_narrow], n_split[i_narrow]))
    print('smallest Omega_max over all orders: %.2f (widest first: %s)' % (omega_max[i_best], i_best == i_wide))
    print('median random order: Omega_max %.2f, shift %.1f pm' % (np.median(omega_max), np.median(shift_max)))
    viol = ok & (shift_max > DELTA * omega_max / (1 - omega_max) + 0.05)
    print('orders violating the bound (21) among Omega<1: %d' % viol.sum())

    FS.apply()
    plt.rcParams.update({'pdf.fonttype': 42})
    fig, ax = plt.subplots(figsize=(3.5, 1.9), layout='constrained')
    rng = np.random.default_rng(1)
    sel = rng.permutation(P)
    single = n_split == 0
    ax.scatter(omega_max[sel][~single[sel]], shift_max[sel][~single[sel]], s=3, c=FS.VERM, lw=0, alpha=.3, rasterized=True, label='split return (%d orders)' % (~single).sum())
    ax.scatter(omega_max[sel][single[sel]], shift_max[sel][single[sel]], s=4, c=FS.BLUE, lw=0, alpha=.9, rasterized=True, label='single peak (%d orders)' % single.sum(), zorder=4)
    om = np.linspace(0.3, 0.97, 200)
    ax.plot(om, DELTA * om / (1 - om), 'k--', lw=.9, label=r'$5\,\Omega_{\max}/(1-\Omega_{\max})$, (21)')
    ax.axvline(1.0, color='0.5', ls=':', lw=.8)
    ax.plot(omega_max[i_wide], shift_max[i_wide], 'o', ms=4.5, mfc='white', mec='k', mew=.8, zorder=6)
    ax.plot(omega_max[i_narrow], shift_max[i_narrow], 's', ms=4.5, mfc='white', mec='k', mew=.8, zorder=6)
    ax.annotate('widest first', (omega_max[i_wide], shift_max[i_wide]), xytext=(6, -2), textcoords='offset points', fontsize=6, va='top')
    ax.annotate('narrowest first', (omega_max[i_narrow], shift_max[i_narrow]), xytext=(-6, -2), textcoords='offset points', fontsize=6, ha='right', va='top')
    ax.text(1.03, 2.3, r'$\Omega_{\max}=1$', fontsize=6, color='0.4')
    ax.set_xlabel(r'Largest $\Omega_k$ of the order')
    ax.set_ylabel('Largest offset of a\nreturn maximum [pm]')
    ax.set_yscale('log')
    ax.set_ylim(2, 400)
    ax.set_xlim(0.6, 3.4)
    ax.legend(fontsize=5.5, loc='lower right', handlelength=1.8, borderaxespad=.3, markerscale=2.5)
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / ('fig_order_permutations.' + ext), dpi=300, bbox_inches='tight', pad_inches=.035)
    print('saved fig_order_permutations')


if __name__ == '__main__':
    main()
