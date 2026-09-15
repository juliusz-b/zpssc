"""s57_estimator_checks.py - the local Gaussian fit against the closed-form rules.

Three checks on Gaussian spectra with the Table III estimator (gauss_fit_peak, +-8 samples
of 20 pm around the largest sample, 30 % threshold):

A(i)  One upstream grating of equal width: fitted shift versus Rule A (whole-line LSQ fit,
      (4/3) sqrt(2/3) R D exp(-D^2/3s^2)) and versus the first-order shift of the return
      maximum, -2 R D exp(-D^2/2s^2). The maxima of the two closed forms differ by exactly 3/2,
      which is the factor 1.5 of the error envelope (34). The local fit lies between them
      (ratio of maxima to Rule A 1.25-1.30 for R = 1-10 %).
A(ii) Co-tuned: n upstream gratings shifted together by delta. Sensitivity of the fitted
      peak versus d lambda_max / d delta = -2nR / (1 - (2n+1)R) of (16). At n = 7, R = 3 %:
      (16) -0.76, exact maximum -0.76, local fit -0.41, whole-line fit -0.26.
B     Leakage at N = 127, W = 200 pm: std of the error at a fixed position grows as sqrt(K)/N,
      the RMS over the band (mean bias included) as K/N and matches c_L K/N with c_L = 27 pm.

python s57_estimator_checks.py
"""
import numpy as np, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import curve_fit
import common as C

SIG = 250.0 / 2.35482          # pm
STEP = 20.0
x = np.arange(-640.0, 640.0 + 1e-9, STEP)   # 65 samples, the Table III grid


def local_fit(y):
    return C.gauss_fit_peak(x, y, win_frac=8.0 / len(x))


def global_fit(y):
    f = lambda t, a, m, s, c: a * np.exp(-0.5 * ((t - m) / s) ** 2) + c
    p, _ = curve_fit(f, x, y, p0=[y.max(), x[np.argmax(y)], SIG, 0.0], maxfev=20000)
    return float(p[1])


print('=== A(i): one upstream grating, equal widths, local fit vs Rule A vs maximum shift')
for R in (0.01, 0.03, 0.10):
    print('R = %.0f%%' % (100 * R))
    print('  Delta/sig   RuleA[pm]  maxshift[pm] (exact)  local[pm]  global[pm]  local/RuleA  local/max')
    for Df in (0.5, 1.0, 1.225, 1.5, 2.0, 2.5):
        D = Df * SIG
        Rk = np.exp(-0.5 * (x / SIG) ** 2)
        Rj = R * np.exp(-0.5 * ((x - D) / SIG) ** 2)
        y = Rk * (1 - Rj) ** 2
        ruleA = -(4 / 3) * np.sqrt(2 / 3) * R * D * np.exp(-D ** 2 / (3 * SIG ** 2))
        mx = -2 * R * D * np.exp(-D ** 2 / (2 * SIG ** 2))
        xf = np.linspace(-300, 300, 600001)
        yf = np.exp(-0.5 * (xf / SIG) ** 2) * (1 - R * np.exp(-0.5 * ((xf - D) / SIG) ** 2)) ** 2
        mxe = xf[np.argmax(yf)]
        loc = local_fit(y); glo = global_fit(y)
        print('  %5.3f   %9.3f   %9.3f (%7.3f) %9.3f  %9.3f     %6.3f     %6.3f' % (Df, ruleA, mx, mxe, loc, glo, loc / ruleA, loc / mx))
    Ds = np.linspace(0.2, 3.0, 57) * SIG
    la, ra = [], []
    for D in Ds:
        Rk = np.exp(-0.5 * (x / SIG) ** 2); Rj = R * np.exp(-0.5 * ((x - D) / SIG) ** 2)
        la.append(abs(local_fit(Rk * (1 - Rj) ** 2)))
        ra.append(abs(-(4 / 3) * np.sqrt(2 / 3) * R * D * np.exp(-D ** 2 / (3 * SIG ** 2))))
    la, ra = np.array(la), np.array(ra)
    print('  max|local| = %.3f pm at Delta/sig = %.2f, max RuleA = %.3f pm, ratio of maxima = %.3f, max pointwise ratio = %.3f'
          % (la.max(), Ds[np.argmax(la)] / SIG, ra.max(), la.max() / ra.max(), (la / ra).max()))

print('\n=== A(ii): co-tuned, n upstream gratings shifted together, sensitivity d(fit)/d(delta)')
for n, R in ((7, 0.01), (7, 0.02), (7, 0.03), (7, 1 / 29), (3, 0.05)):
    sens_th = -2 * n * R / (1 - (2 * n + 1) * R)
    out = []
    for est in (local_fit, global_fit):
        vals = []
        for d in (-4.0, 4.0):
            Rk = np.exp(-0.5 * (x / SIG) ** 2)
            Rj = R * np.exp(-0.5 * ((x - d) / SIG) ** 2)
            vals.append(est(Rk * (1 - Rj) ** (2 * n)))
        out.append((vals[1] - vals[0]) / 8.0)
    vals = []
    for d in (-4.0, 4.0):
        xf = np.linspace(-100, 100, 400001)
        yf = np.exp(-0.5 * (xf / SIG) ** 2) * (1 - R * np.exp(-0.5 * ((xf - d) / SIG) ** 2)) ** (2 * n)
        vals.append(xf[np.argmax(yf)])
    mx = (vals[1] - vals[0]) / 8.0
    print('  n=%d R=%.4f  (16): %.3f   exact max: %.3f   local fit: %.3f   global fit: %.3f' % (n, R, sens_th, mx, out[0], out[1]))

print('\n=== B: leakage error scaling, N = 127, W = 200 pm, sigma = 106 pm, local fit')
N = 127
_ms = 1.0 - 2.0 * C._mls01(7)
AC = C.periodic_xcorr(_ms, _ms)
W = 200.0
rng = np.random.default_rng(0)
print('   K   std|nu_k=0   std|nu_k=150   RMS over band   27*K/N')
for K in (4, 8, 16, 32, 64):
    def trial(nu_fixed=None, ntr=300):
        errs = []
        for t in range(ntr):
            bins = np.sort(rng.choice(np.arange(1, N), size=K, replace=False))
            nub = rng.uniform(-W, W, size=K)
            if nu_fixed is not None:
                nub[0] = nu_fixed
            A = np.exp(-0.5 * ((x[None, :] - nub[:, None]) / SIG) ** 2)
            Wl = AC[(bins[:, None] - bins[None, :]) % N]
            np.fill_diagonal(Wl, 0.0)
            S = A + Wl @ A
            if nu_fixed is not None:
                errs.append(local_fit(S[0]) - nub[0])
            else:
                errs.extend([local_fit(S[k]) - nub[k] for k in range(K)])
        return np.array(errs)
    e0 = trial(0.0, 200); e1 = trial(150.0, 200); eb = trial(None, max(20, 400 // K))
    print('  %3d   %7.2f      %7.2f        %7.2f        %6.2f' % (K, e0.std(), e1.std(), np.sqrt(np.mean(eb ** 2)), 27.0 * K / N))
