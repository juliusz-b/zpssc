"""s42_taper.py - does a reflectivity profile along the fiber help?

Shadowing on grating k comes only from the gratings in front of it, so the
error depends on where the strong gratings sit, not only on how strong they
are on average. Four profiles with the same mean reflectivity are compared:

  uniform     every grating at the mean
  rising      weak at the interrogator end, strong at the far end
  falling     strong at the interrogator end, weak at the far end
  equalized   the classical TDM profile, R_k chosen so that every direct
              return has the same power at the detector

Shadowing-only RMS error over the same random detuning draws for every
profile, with and without the sequential correction of (13), and the weakest
direct return in dB relative to the mean reflectivity.
"""
import numpy as np
import common as C

K = 16
FWHM = 250.0
SIG = FWHM / 2.35482
NDRAW = 8
SEED = 3


def line(nu, nu0):
    return np.exp(-0.5 * ((nu - nu0) / SIG) ** 2)


def equalized(mean):
    """R_k / prod_{j<k}(1-R_j)^2 constant, scaled to the requested mean."""
    def build(r1):
        r = np.empty(K)
        t = 1.0
        for k in range(K):
            r[k] = min(r1 / t, 0.95)
            t *= (1.0 - r[k]) ** 2
        return r
    lo, hi = 1e-5, mean
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if build(mid).mean() < mean:
            lo = mid
        else:
            hi = mid
    return build(0.5 * (lo + hi))


def profiles(mean):
    lin = np.linspace(0.2, 1.8, K) * mean
    return {'uniform': np.full(K, mean), 'rising': lin,
            'falling': lin[::-1], 'equalized': equalized(mean)}


def run(R, peel):
    rng = np.random.default_rng(SEED)
    nu = np.linspace(-520.0, 520.0, 1041)
    errs = []
    for _ in range(NDRAW):
        cen = rng.uniform(-200.0, 200.0, K)
        shapes = np.array([line(nu, c) for c in cen])
        trans = np.ones((K, nu.size))
        for k in range(1, K):
            trans[k] = trans[k - 1] * (1.0 - R[k - 1] * shapes[k - 1]) ** 2
        A = R[:, None] * shapes * trans
        if peel:
            est = np.ones_like(nu)
            corr = np.empty_like(A)
            for k in range(K):
                corr[k] = A[k] / np.maximum(est, 0.05)
                est = est * (1.0 - np.clip(corr[k], 0.0, 0.99)) ** 2
            A = corr
        for k in range(K):
            errs.append(C.gauss_fit_peak(nu, A[k]) - cen[k])
    return np.sqrt(np.mean(np.square(errs)))


def weakest_db(R):
    """Direct return of the last grating when everything is co-tuned."""
    t = np.cumprod(np.concatenate([[1.0], (1.0 - R[:-1]) ** 2]))
    return 10.0 * np.log10((R * t).min() / R.mean())


for mean in (0.10, 0.05, 0.02):
    print('--- K = %d, mean reflectivity %.0f%% ---' % (K, 100 * mean))
    print('   %-10s %10s %12s %14s  profile' % (
        'profile', 'raw [pm]', 'peeled [pm]', 'weakest [dB]'))
    for name, R in profiles(mean).items():
        print('   %-10s %10.2f %12.2f %14.1f  %s' % (
            name, run(R, False), run(R, True), weakest_db(R),
            ' '.join('%.1f' % (100 * r) for r in R[:4])
            + ' ... ' + ' '.join('%.1f' % (100 * r) for r in R[-3:])))
