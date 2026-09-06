"""s61_crb_sampling.py - how many sweep steps does a Bragg-wavelength estimate
need? Cramer-Rao bound for the centre of a Gaussian line.

Model: y_i = a exp(-(x_i - mu)^2 / 2 s^2) + b + n_i, white noise sigma_n,
unknowns (mu, a, s, b) or a subset. Fisher matrix J = G' G / sigma_n^2 with
G the Jacobian at the sample points. The bound on mu is [J^-1]_mumu.

Two questions:
  1. For M equally spaced steps over the band, how does the bound fall with M
     (it goes as 1/M, dense sampling only averages noise)?
  2. Where should a few samples sit? For mu alone the Fisher weight of a
     sample is (x/s^2)^2 exp(-x^2/s^2), maximal at x = +-s (the inflection
     points, the steepest slope). With a, s, b unknown as well the optimum
     moves, so we search it numerically for 3, 4, 5 and 6 points.

Then the design consequence: with the refresh rate of (23) proportional to
1/M, the same precision at a fraction of the steps, or the same time spent
on more samples per step.

Output: out/s61_crb_sampling.txt
"""
import os
import numpy as np
from scipy.optimize import minimize

s = 1.0                      # line width parameter, everything in units of sigma
a, b = 1.0, 0.0


def jac(x, mu=0.0):
    e = np.exp(-0.5 * ((x - mu) / s) ** 2)
    return np.column_stack([a * e * (x - mu) / s ** 2,     # d/dmu
                            e,                              # d/da
                            a * e * (x - mu) ** 2 / s ** 3,  # d/ds
                            np.ones_like(x)])               # d/db


def crb_mu(x, unknowns):
    G = jac(np.asarray(x, float))[:, unknowns]
    J = G.T @ G
    try:
        return float(np.linalg.inv(J)[0, 0])       # in units of (sigma_n/a)^2 * s^2
    except np.linalg.LinAlgError:
        return np.inf


lines = []
say = lambda t='': (print(t), lines.append(t))
say('Cramer-Rao bound on the fitted centre, in units of (sigma_n/a) x sigma_line')
say('unknowns: (mu, a, s, b) unless stated')
U4 = [0, 1, 2, 3]; U1 = [0]
say()
say('--- equally spaced steps over +-3 sigma (the 64-step sweep of Table III covers +-6 sigma of a 250-pm line) ---')
say('%6s %14s %14s' % ('M', 'sqrt(CRB) 4 unk', 'sqrt(CRB) mu only'))
for M in (4, 5, 6, 8, 12, 16, 32, 64):
    x = np.linspace(-3, 3, M)
    say('%6d %14.3f %14.3f' % (M, np.sqrt(crb_mu(x, U4)), np.sqrt(crb_mu(x, U1))))
say()
say('--- best placement of a few samples (numerical search, symmetric start) ---')
for n in (3, 4, 5, 6):
    best = None
    for trial in range(40):
        x0 = np.sort(np.random.default_rng(trial).uniform(-2.5, 2.5, n))
        r = minimize(lambda x: crb_mu(x, U4) if len(set(np.round(x, 3))) == n else 1e9, x0, method='Nelder-Mead',
                     options=dict(xatol=1e-4, fatol=1e-8, maxiter=4000))
        if np.isfinite(r.fun) and (best is None or r.fun < best.fun):
            best = r
    say('%d samples: sqrt(CRB) = %.3f at x/sigma = %s' % (n, np.sqrt(best.fun), np.round(np.sort(best.x), 2)))
say()
say('--- mu only (a, s, b known): two samples at +-1 sigma give sqrt(CRB) = %.3f, the equally spaced 64-step sweep gives %.3f' %
    (np.sqrt(crb_mu([-1.0, 1.0], U1)), np.sqrt(crb_mu(np.linspace(-3, 3, 64), U1))))
say()
say('--- reading: the 64-step sweep spends most samples where the line carries no information about its centre.')
say('    Four samples placed by the search reach the precision of about %d equally spaced steps,' %
    int(round((np.sqrt(crb_mu(np.linspace(-3, 3, 64), U4)) / np.sqrt(crb_mu(minimize(lambda x: crb_mu(x, U4), [-1.5, -0.5, 0.5, 1.5], method="Nelder-Mead").x, U4))) ** 2 * 64)))
say('    so a tracking sweep of a few steps per line raises the refresh rate of (23) by that factor at equal precision,')
say('    or keeps M and averages that many more code periods per step.')
os.makedirs('out', exist_ok=True)
open('out/s61_crb_sampling.txt', 'w').write('\n'.join(lines) + '\n')
