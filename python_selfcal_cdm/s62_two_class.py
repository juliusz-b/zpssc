"""s62 - two-class despreading: the code's transitions as a second wavelength probe.

An on-chip that follows a zero carries the transient chirp overshoot, an on-chip
that follows a one sits at the adiabatic level. The code therefore splits into
two sub-codes, c_v (after zero) and c_w (after one), and every grating returns
two spectra, one per class, separated on the wavelength axis by the class-mean
chirp difference. Least squares against the shifted sub-codes recovers both
spectra per grating per sweep step. The shift between their fitted centres
measures the pattern-dependent chirp step in situ, and the bit-history echoes
of the ordinary correlator (lags 0, 1 and the shift-and-add lag s) disappear.

Time-domain model as in s55: unipolar m-sequence, SPC samples per chip,
xi(t) = A_ad on every on-chip plus A_T exp(-t/tau) after each 0->1 edge,
Gaussian line, 64-step sweep, chip-averaged record, no noise.

Figure (single column, two panels):
  (a) the two recovered spectra of one grating, their difference, and an inset
      of recovered versus true chirp step for five chirp waveforms;
  (b) centre bias of the grating under test versus the lag of a neighbour,
      ordinary correlator versus two-class decoder.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import max_len_seq
import common as C
import figstyle as FS
import warnings
warnings.filterwarnings('ignore')

NB, SPC = 7, 32
N = 2 ** NB - 1
FWHM = 250.0
SIG = FWHM / 2.35482
u = max_len_seq(NB)[0].astype(float)
b = 2 * u - 1
v = u * (1 - np.roll(u, 1))
w = u * np.roll(u, 1)
drive = np.repeat(u, SPC)
L = drive.size
LAM = np.linspace(-650, 650, 64)
os.makedirs('out', exist_ok=True)
os.makedirs('figs', exist_ok=True)
lines = []


def say(s=''):
    print(s)
    lines.append(s)


def xi_wave(A_ad, A_T, tau_chip):
    x = np.full(L, A_ad, float)
    for e in np.where(np.diff(np.r_[u[-1], u]) > 0)[0]:
        i0 = e * SPC
        x[i0:] += A_T * np.exp(-np.arange(L - i0) / (tau_chip * SPC))
    return x


def record(gratings, xi, lam):
    y = np.zeros(L)
    for dd, dl, amp in gratings:
        sh = int(round(dd * SPC))
        y += amp * np.roll(drive, sh) * np.exp(-0.5 * ((lam - dl + np.roll(xi, sh)) / SIG) ** 2)
    yc = y.reshape(N, SPC).mean(axis=1)
    return yc - yc.mean()


def cols(delays, classes):
    cs = []
    for dd in delays:
        for cl in classes:
            c = np.roll(cl, int(round(dd)))
            cs.append(c - c.mean())
    return np.array(cs).T


def class_mean(xi, cl):
    m = np.repeat(cl, SPC) > 0
    return xi[m].mean()


def spectra(gratings, xi, delays, classes):
    P = np.linalg.pinv(cols(delays, classes))
    return np.array([P @ record(gratings, xi, lam) for lam in LAM])


say('s62 two-class despreading, N=%d, SPC=%d, FWHM=%.0f pm' % (N, SPC, FWHM))

# --- panel (a): one grating, two spectra, identity of the step ------------------
WAVES = ((20, 30, 0.1), (20, 30, 0.3), (20, 30, 1.0), (20, 80, 0.3), (0, 80, 0.3))
steps_true, steps_fit = [], []
say('\n(a) single grating: true class-mean step vs shift of the fitted centres [pm]')
for A_ad, A_T, tau in WAVES:
    xi = xi_wave(A_ad, A_T, tau)
    A = spectra([(0, 0.0, 1.0)], xi, [0], [v, w])
    cv, cw = C.gauss_fit_peak(LAM, A[:, 0]), C.gauss_fit_peak(LAM, A[:, 1])
    t = class_mean(xi, v) - class_mean(xi, w)
    steps_true.append(t)
    steps_fit.append(cw - cv)
    say('   A_ad=%3d A_T=%3d tau=%.1f chip: true %6.2f  recovered %6.2f' % (A_ad, A_T, tau, t, cw - cv))
xi_show = xi_wave(20, 30, 0.3)
A_show = spectra([(0, 0.0, 1.0)], xi_show, [0], [v, w])
Av, Aw = A_show[:, 0], A_show[:, 1]
peak = max(Av.max(), Aw.max())
Av, Aw = Av / peak, Aw / peak
step_show = steps_fit[1]

# --- panel (b): neighbour lag sweep ---------------------------------------------
say('\n(b) grating under test at delay 0, neighbour d chips closer, A_T=30 pm tau=0.3 chip, detuning 0: centre bias [pm]')
xi = xi_wave(20, 30, 0.3)
xi0 = np.full(L, class_mean(xi, u))
LAGS = np.arange(1, 13)
bias_mf, bias_tc = [], []
for d in LAGS:
    g = [(0, 0.0, 1.0), (-int(d), 0.0, 1.0)]
    Su = spectra(g, xi, [0, -int(d)], [u])[:, 0]
    Sr = spectra(g, xi0, [0, -int(d)], [u])[:, 0]
    A = spectra(g, xi, [0, -int(d)], [v, w])
    Stc = 0.5 * (A[:, 0] + A[:, 1])
    cu, cr, ct = (C.gauss_fit_peak(LAM, s) for s in (Su, Sr, Stc))
    bias_mf.append(cu - cr)
    bias_tc.append(ct - cr)
    say('   d=%2d  correlator %6.2f  two-class %6.2f' % (d, cu - cr, ct - cr))
bias_mf, bias_tc = np.array(bias_mf), np.array(bias_tc)
s_lag = [l for l in range(2, N) if abs((np.roll(b, l) * v).sum() / (u * b).sum()) > 0.1][0]
say('   shift-and-add lag of this sequence: s=%d' % s_lag)

# --- conditioning versus lag -------------------------------------------------
say('\n(c) smallest singular value of [D_v D_w], normalized, two gratings at lag d')
for d in (1, 2, 6, 7, 8):
    sv = np.linalg.svd(cols([0, d], [v, w]), compute_uv=False)
    say('   d=%d: %.3f' % (d, sv[-1] / sv[0]))

# --- figure ---------------------------------------------------------------------
FS.apply()
fig, ax = plt.subplots(2, 1, figsize=(3.45, 4.6))
a = ax[0]
a.plot(LAM, Av, color=FS.BLUE, label='After a zero, $A^{v}$')
a.plot(LAM, Aw, color=FS.VERM, ls='--', label='After a one, $A^{w}$')
a.plot(LAM, 10 * (Aw - Av), color=FS.GREY, lw=1.2, label='$10\\,(A^{w}-A^{v})$')
a.axhline(0, color='k', lw=0.6)
a.set_xlim(-650, 650)
a.set_ylim(-0.65, 1.15)
a.set_xlabel('Wavelength offset (pm)')
a.set_ylabel('Normalized amplitude')
a.legend(loc='lower left', fontsize=6.5)
a.text(0.97, 0.05, 'Fitted-center shift %.1f pm' % step_show, transform=a.transAxes, ha='right', va='bottom', fontsize=7)
FS.letter(a, 'a')
ins = a.inset_axes([0.66, 0.56, 0.31, 0.38])
mx = max(steps_true) * 1.15
ins.plot([0, mx], [0, mx], color=FS.LGREY, lw=0.8)
ins.plot(steps_true, steps_fit, 'o', color=FS.GREEN, ms=3.5)
ins.set_xlim(0, mx)
ins.set_ylim(0, mx)
ins.set_xlabel('True step (pm)', fontsize=6, labelpad=1)
ins.set_ylabel('Recovered (pm)', fontsize=6, labelpad=1)
ins.tick_params(labelsize=5.5, length=2)
ins.set_xticks([0, 10, 20])
ins.set_yticks([0, 10, 20])

bx = ax[1]
bx.plot(LAGS, bias_mf, 'o-', color=FS.ORANGE, label='Correlator')
bx.plot(LAGS, bias_tc, 's-', color=FS.BLUE, label='Two-class decoder')
bx.axhline(0, color='k', lw=0.6)
bx.set_xlim(0.5, 12.5)
bx.set_ylim(-5.5, 5.5)
bx.set_xticks(LAGS)
bx.set_xlabel('Neighbor lag (chips)')
bx.set_ylabel('Center bias (pm)')
bx.legend(loc='upper right', fontsize=6.5)
bx.annotate('lag 1', (1, bias_mf[0]), xytext=(1.6, bias_mf[0] - 0.25), fontsize=6.5)
bx.annotate('lag $s$ = %d' % s_lag, (s_lag, bias_mf[s_lag - 1]), xytext=(s_lag + 0.6, bias_mf[s_lag - 1] + 0.3), fontsize=6.5)
FS.letter(bx, 'b')
fig.tight_layout(h_pad=1.0)
for ext in ('pdf', 'png'):
    fig.savefig('figs/fig_s62_twoclass.' + ext, dpi=300)
open('out/s62_two_class.txt', 'w').write('\n'.join(lines) + '\n')
say('\nfigure: figs/fig_s62_twoclass.pdf')
