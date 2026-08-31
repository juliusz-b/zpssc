"""s44_additivity.py - do the pairwise Law A shifts add up?

Law A prices one upstream neighbour. A real array has many, and the paper
claims that to first order in R the pairwise shifts add, so the sum over
neighbours tracks the full serial product of (8). This checks that claim:
random detunings, the fitted peak of R_k * prod (1 - R_j)^2 against the sum
of the closed-form pairwise terms, for growing n and both reflectivities.

Expected outcome, quoted in the paper: at R = 1 percent the sum tracks the
full model to about a tenth up to one hundred upstream gratings. At
R = 10 percent it stops being quantitative a few gratings deep, but the
attenuation kills the measurement there first anyway: one hundred 10
percent gratings transmit 0.9^200, 92 dB down.
"""
import warnings

import numpy as np
from scipy.optimize import curve_fit

warnings.filterwarnings('ignore')

SIG = 106.2
CA = 4.0 / 3.0 * np.sqrt(2.0 / 3.0)
NDRAW = 6

nu = np.linspace(-5.0 * SIG, 5.0 * SIG, 1800)
wanted = np.exp(-0.5 * (nu / SIG) ** 2)


def gauss(x, a, m, s, b):
    return b + a * np.exp(-0.5 * ((x - m) / s) ** 2)


def fitpeak(y):
    p, _ = curve_fit(gauss, nu, y, p0=[0.9, 0.0, SIG, 0.0], maxfev=20000)
    return p[1]


rng = np.random.default_rng(1)
print('%6s %5s %28s %28s %10s' % ('R', 'n', 'full serial model [pm]',
                                  'sum of Law A terms [pm]',
                                  'co-tuned'))
for R, n in ((0.10, 5), (0.10, 10), (0.01, 10), (0.01, 50), (0.01, 100)):
    full_v, add_v = [], []
    for _ in range(NDRAW):
        dets = rng.uniform(-200.0, 200.0, n)
        prod = np.ones_like(nu)
        for d in dets:
            prod *= (1.0 - R * np.exp(-0.5 * ((nu - d) / SIG) ** 2)) ** 2
        full_v.append(fitpeak(wanted * prod))
        add_v.append(sum(-CA * R * d * np.exp(-d ** 2 / (3.0 * SIG ** 2))
                         for d in dets))
    att = 10.0 * np.log10((1.0 - R) ** (2 * n))
    print('%5.0f%% %5d %28s %28s %8.1f dB'
          % (100 * R, n,
             ' '.join('%+7.2f' % v for v in full_v[:3]),
             ' '.join('%+7.2f' % v for v in add_v[:3]), att))

print()
print('transmission after 100 gratings, co-tuned worst case:')
for R in (0.10, 0.01):
    t = (1.0 - R) ** 200
    print('  R = %2.0f%%: %g  (%.1f dB)' % (100 * R, t, 10 * np.log10(t)))
