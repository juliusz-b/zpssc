"""s_fig1_principle.py - Fig. 1 of the paper: stepped-sweep CDM interrogation of four gratings.

(a) optical path with the gratings at z_k and the measured VCSEL tuning curve, (b) the four lines sharing one band
and the samples A_k(lambda_m) taken at one step, (c) HCG voltage steps, bias + code current and the optical output,
(d) delayed echoes at the photodiode and their sum P_m(t), (e) correlation X_m(tau) at one step with the peaks at
tau_1..tau_4, (f) the spectrum S_2(lambda_m) assembled from the M steps and its fitted center.
Everything is illustrative (K = 4, R = 10 %, N = 127, M = 25). The tuning curve is the quartic fit of the measured
HCG-VCSEL characteristic used by the decoder of the paper.

Output: figs/fig_s_fig1_principle.{pdf,png}
"""
import sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyBboxPatch
import common as C


class D:
    """Measured HCG-VCSEL tuning curve: quartic fit lambda(V) [nm] for V in 0..14 V, the same numbers as in the decoder."""
    TUNE_P4 = np.array([-7.14581756e-05, 1.94718096e-03, -4.51295135e-02, -9.75634596e-02, 1.56974137e+03])
    TUNE_V = np.linspace(0.0, 14.0, 2801)
    TUNE_L = np.polyval(TUNE_P4, TUNE_V)


import figstyle as FS
FS.apply()
plt.rcParams.update({'pdf.fonttype': 42})
from pathlib import Path
OUT = Path(__file__).resolve().parent / "figs"

COLS = [FS.VERM, FS.ORANGE, FS.GREEN, FS.PURPLE]
K, R, N, M = 4, 0.10, 127, 25
SIG = 250.0 / 2.35482                        # pm
DET = np.array([-90.0, 30.0, -30.0, 90.0])   # Bragg offsets of FBG1..4 from the band center [pm]
Z = np.array([4.0, 11.5, 21.0, 35.0])        # positions [m]
DZ = 4.09                                    # one chip at 25 Mchip/s [m]
TAU = Z / DZ                                 # delays in chips
LAM_M = 0.0                                  # the wavelength step shown in (b), (d), (e)
A = R * np.exp(-0.5 * ((LAM_M - DET) / SIG) ** 2)   # A_k(lambda_m)
code = C._mls01(7).astype(float)             # N = 127 chips, 0/1

fig = plt.figure(figsize=(7.1, 3.45), layout='constrained')
gs = fig.add_gridspec(2, 3, width_ratios=[1.55, 1.0, 1.0], height_ratios=[1.3, 1.0])
axa = fig.add_subplot(gs[0, 0]); axb = fig.add_subplot(gs[0, 1]); axc = fig.add_subplot(gs[0, 2])
axd = fig.add_subplot(gs[1, 0]); axe = fig.add_subplot(gs[1, 1]); axf = fig.add_subplot(gs[1, 2])

# ---------------- (a) optical path ----------------
ax = axa; ax.axis('off'); ax.set_xlim(0, 10); ax.set_ylim(0.3, 4.95)
def block(x, y, w, h, t, fc='#F2F2F2'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.02,rounding_size=0.12', fc=fc, ec='0.25', lw=0.8))
    ax.text(x + w / 2, y + h / 2, t, ha='center', va='center', fontsize=7)
YF = 3.55                                    # fiber line
block(0.1, YF - 0.45, 1.6, 0.9, 'VCSEL')
block(0.1, 0.55, 1.6, 0.9, 'PD')
ax.add_patch(Circle((2.7, YF), 0.42, fc='white', ec='0.25', lw=0.8)); ax.text(2.7, YF, 'circ.', ha='center', va='center', fontsize=6.2)
ax.plot([1.7, 2.28], [YF, YF], color='0.2', lw=1.3)
ax.plot([2.7, 2.7, 1.95], [YF - 0.42, 1.0, 1.0], color='0.2', lw=1.3)
ax.annotate('', xy=(1.72, 1.0), xytext=(2.2, 1.0), arrowprops=dict(arrowstyle='-|>', color='0.2', lw=1.1, mutation_scale=8))
ax.plot([3.12, 9.9], [YF, YF], color='0.2', lw=1.8)
x0, scale = 3.5, (9.65 - 3.5) / Z.max()
for k in range(K):
    xk = x0 + Z[k] * scale
    ax.add_patch(Rectangle((xk - 0.25, YF - 0.42), 0.5, 0.84, fc=COLS[k], ec='none'))
    ax.text(xk, YF - 0.5, 'FBG%d' % (k + 1), ha='center', va='top', fontsize=6, color=COLS[k])
    y = YF + 0.6 + 0.25 * k
    ax.annotate('', xy=(xk, y), xytext=(3.12, y), arrowprops=dict(arrowstyle='<->', color=COLS[k], lw=0.8, mutation_scale=5, shrinkA=0, shrinkB=0))
    ax.text(xk + 0.14, y, '$z_%d$' % (k + 1), ha='left', va='center', fontsize=6.5, color=COLS[k])
ax.text(3.3, 2.58, 'round trip $\\tau_k=2n_gz_k/c$', fontsize=6.5, ha='left', va='center', color='0.25')
ax.text(3.3, 2.26, 'FBG$k$ is read through the gratings before it', fontsize=6, ha='left', va='center', color='0.4')
# inset: measured VCSEL tuning curve, the four steps of (c) marked on it
axi = ax.inset_axes([0.60, 0.02, 0.39, 0.34])
axi.plot(D.TUNE_V, D.TUNE_L, color='0.35', lw=1.0)
V4 = np.array([3.0, 5.0, 7.0, 9.0])
axi.plot(V4, np.polyval(D.TUNE_P4, V4), 'o', color=FS.PURPLE, ms=2.8, mec='white', mew=0.4, zorder=5)
for i, v in enumerate(V4):
    axi.text(v - 0.4, np.polyval(D.TUNE_P4, v) - 0.6, '$\\lambda_%d$' % (i + 1), fontsize=5.8, color=FS.PURPLE, ha='right', va='top')
axi.set_xlim(0, 14); axi.set_ylim(1561.5, 1571.5)
axi.set_xticks([0, 7, 14]); axi.set_yticks([1562, 1566, 1570])
axi.tick_params(labelsize=5.5, length=2, pad=1)
axi.set_xlabel('$V_{\\mathrm{HCG}}$ [V]', fontsize=5.5, labelpad=0.5); axi.set_ylabel('$\\lambda$ [nm]', fontsize=5.5, labelpad=0.5)
axi.text(0.96, 0.92, 'VCSEL tuning', transform=axi.transAxes, fontsize=5.2, ha='right', va='top', color='0.35')
FS.letter(ax, 'a')

# ---------------- (b) four lines in one band ----------------
ax = axb
lam = np.linspace(-450, 450, 600)
for k in range(K):
    ax.plot(lam, np.exp(-0.5 * ((lam - DET[k]) / SIG) ** 2), color=COLS[k], lw=1.2, label='FBG%d' % (k + 1))
ax.axvline(LAM_M, color='0.4', lw=0.8, ls=(0, (2, 2)))
for k in range(K):
    ax.plot(LAM_M, A[k] / R, 'o', color=COLS[k], ms=3.2, mec='white', mew=0.5, zorder=5)
ax.text(LAM_M + 15, 1.16, '$\\lambda_m$', fontsize=6.5, ha='left', va='center')
ax.set_xlim(-450, 450); ax.set_ylim(0, 1.45); ax.set_yticks([0, 0.5, 1.0])
ax.set_xlabel('Wavelength offset [pm]'); ax.set_ylabel('$R_k(\\lambda)/R$')
ax.legend(fontsize=5.5, loc='upper center', ncol=4, handlelength=1.0, columnspacing=0.5, handletextpad=0.3, borderaxespad=0.2, frameon=False)
FS.letter(ax, 'b')

# ---------------- (c) drive ----------------
ax = axc
t = np.linspace(0, 4, 2000)
step = np.floor(t).clip(0, 3)
chips = code[(np.floor(t * 12) % 12).astype(int)]     # 12 chips per step shown
Y_V, Y_I, Y_P = 2.75, 1.5, 0.25
ax.plot(t, Y_V + 0.17 * step, color=FS.PURPLE, lw=1.2)
for s in range(4):
    ax.text(s + 0.5, Y_V + 0.17 * s + 0.09, '$\\lambda_%d$' % (s + 1), ha='center', va='bottom', fontsize=6, color=FS.PURPLE)
ax.plot(t, Y_I + 0.5 * chips, color=FS.BLUE, lw=0.9)
ax.plot(t, Y_P + 0.5 * chips * (0.85 + 0.05 * step), color='0.2', lw=0.9)
for s in range(1, 4):
    ax.axvline(s, color='0.75', lw=0.5, ls=(0, (2, 2)))
ax.set_yticks([Y_P + 0.2, Y_I + 0.2, Y_V + 0.25])
ax.set_yticklabels(['$P_{\\mathrm{out}}$', '$I$', '$V_{\\mathrm{HCG}}$'])
ax.tick_params(axis='y', length=0)
ax.set_xlim(0, 4); ax.set_ylim(0, 3.75)
ax.set_xticks([0.5, 1.5, 2.5, 3.5]); ax.set_xticklabels(['step 1', '2', '3', '4'])
ax.set_xlabel('Time')
ax.text(3.95, Y_I + 0.62, 'bias $+$ code $c(t)$', fontsize=5.4, ha='right', va='bottom', color=FS.BLUE)
ax.text(3.95, Y_P + 0.62, 'light at $\\lambda_m$ carries $c(t)$', fontsize=5.4, ha='right', va='bottom', color='0.2')
FS.letter(ax, 'c')

# ---------------- (d) echoes and their sum ----------------
ax = axd
tt = np.linspace(0, 16, 3200)
def echo(k):
    idx = np.floor(tt - TAU[k]).astype(int)
    out = np.where(idx >= 0, code[np.clip(idx, 0, N - 1) % N], 0.0)
    return A[k] / R * out
tot = np.zeros_like(tt)
for k in range(K):
    e = echo(k); tot += e
    ax.plot(tt, 3.9 - 0.85 * k + 0.55 * e, color=COLS[k], lw=0.8)
    ax.text(16.3, 3.9 - 0.85 * k + 0.25, 'FBG%d: $A_%d(\\lambda_m)\\,c(t-\\tau_%d)$' % (k + 1, k + 1, k + 1), fontsize=5.5, va='center', ha='left', color=COLS[k])
ax.plot(tt, 0.05 + 0.32 * tot, color='0.15', lw=0.9)
ax.text(16.3, 0.35, 'sum: $P_m(t)$', fontsize=5.5, va='center', ha='left', color='0.15')
yb = 3.9 - 0.85 * 1 - 0.18
ax.annotate('', xy=(TAU[1], yb), xytext=(0, yb), arrowprops=dict(arrowstyle='-|>', color=COLS[1], lw=0.7, mutation_scale=6, shrinkA=0, shrinkB=0))
ax.text(TAU[1] / 2, yb - 0.05, '$\\tau_2$', fontsize=6, ha='center', va='top', color=COLS[1])
ax.set_xlim(0, 16); ax.set_ylim(-0.2, 4.7); ax.set_yticks([]); ax.tick_params(axis='y', length=0)
ax.set_xlabel('Time [chips]'); ax.set_ylabel('Optical power')
FS.letter(ax, 'd')

# ---------------- (e) correlation at one wavelength step ----------------
ax = axe
tau = np.linspace(0, 12, 1200)
rng = np.random.default_rng(3)
X = -np.sum(A) / N + 0.0015 * rng.normal(size=tau.size)
for k in range(K):
    X = X + A[k] * np.clip(1 - np.abs(tau - TAU[k]), 0, None)
ax.plot(tau, X, color='0.25', lw=0.9)
for k in range(K):
    ax.plot(TAU[k], A[k], 'o', color=COLS[k], ms=3.4, mec='white', mew=0.5, zorder=5)
    ax.text(TAU[k], -0.012, '$\\tau_%d$' % (k + 1), ha='center', va='top', fontsize=6, color=COLS[k])
ax.text(TAU[1] + 0.35, A[1] + 0.003, '$S_2(\\lambda_m)$', fontsize=6, ha='left', va='bottom', color=COLS[1])
ax.set_xlim(0, 12); ax.set_ylim(-0.03, 0.125); ax.set_yticks([0, 0.05, 0.1])
ax.set_xlabel('Delay $\\tau$ [chips]'); ax.set_ylabel('$X_m(\\tau)$')
FS.letter(ax, 'e')

# ---------------- (f) spectrum of grating 2 ----------------
ax = axf
lm = np.linspace(-240, 240, M)
S2 = R * np.exp(-0.5 * ((lm - DET[1]) / SIG) ** 2) + 0.0015 * rng.normal(size=M)
lf = np.linspace(-240, 240, 400)
ax.plot(lf, R * np.exp(-0.5 * ((lf - DET[1]) / SIG) ** 2), color=COLS[1], lw=1.0, alpha=0.6)
ax.plot(lm, S2, 'o', color=COLS[1], ms=2.6, mec='white', mew=0.4)
i_m = np.argmin(np.abs(lm - LAM_M))
ax.plot(LAM_M, S2[i_m], 'o', color=COLS[1], ms=4.2, mec='0.2', mew=0.8, zorder=6)
ax.axvline(DET[1], color='0.4', lw=0.8, ls=(0, (2, 2)))
ax.text(DET[1] + 8, 0.02, '$\\hat\\lambda_{B,2}$', fontsize=6.5, ha='left', va='center')
ax.text(LAM_M - 10, S2[i_m] + 0.004, 'step $m$', fontsize=5.5, ha='right', va='bottom', color='0.3')
ax.set_xlim(-240, 240); ax.set_ylim(-0.005, 0.125); ax.set_yticks([0, 0.05, 0.1])
ax.set_xlabel('$\\lambda_m$ offset [pm]'); ax.set_ylabel('$S_2(\\lambda_m)$')
FS.letter(ax, 'f')

OUT.mkdir(exist_ok=True)
fig.savefig(OUT / 'fig_s_fig1_principle.png', dpi=300); fig.savefig(OUT / 'fig_s_fig1_principle.pdf')
print('saved figs/fig_s_fig1_principle.pdf and .png')
print('ok')
