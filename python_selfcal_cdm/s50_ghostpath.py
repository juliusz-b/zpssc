"""s50_ghostpath.py - every arrival at the photodiode for three gratings.

Replaces the TikZ drawing fig_ghost.tex with a matplotlib figure in the house
style (figstyle), so that fonts, colours and panel letters match the other
figures of the paper.

  (a) Space-time diagram. Horizontal axis: time since the launch, in units of
      the delay. Vertical axis: position along the fibre, with the fibre
      itself, its three gratings 1, 2, 3 and the photodiode drawn on the
      left. The launched code climbs the diagonal; every reflection descends
      back to z = 0 and lands on the time axis at its delay. The three direct
      returns are blue, and the third-order path (3,1,2) is highlighted in
      orange. Other paths appear in the arrival strip of panel (b).
  (c) Spectrum of the last grating of an eight-grating uniform array at
      R = 10 % (the spectral model of s41): shadowing lowers and shifts the
      measured line, the sequential correction restores it, but the ghosts
      under this grating stay and are scaled up together with the line.
  (b) The delay profile of the same three gratings in the time-domain model:
      z = 8, 19.2, 40 m (the ratios of (a)), R = 10 %, unipolar m-sequence
      N = 127 at 25 Mchip/s (Table III), every path up to the third order,
      record mean removed, correlation with the bipolar replica, constant
      offset of the periodic correlation subtracted, linear scale. Above the
      curve, a separate strip marks three direct returns with squares and
      five third-order paths with circles. Stacking indicates paths at the
      same delay, not optical power. An inset magnifies the late echoes.
      The spacing is non-uniform on purpose so that the
      ghosts do not collide, except the pair (3,1,2)/(2,1,3), which share one
      delay for any spacing.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from matplotlib.colors import to_rgb
import common as C
import figstyle as FS

FS.apply()

# ---------------------------------------------------------------------------
# geometry: time in delay units, position z = tau / 2 (c / n_g = 1)
# ---------------------------------------------------------------------------
TB, TC, TA = 1.0, 2.4, 5.0
T = {'a': TA, 'b': TB, 'c': TC}
COL = {'a': FS.VERM, 'b': FS.ORANGE, 'c': FS.GREEN}
GREY = '0.45'
NAME = {'b': '1', 'c': '2', 'a': '3'}   # gratings numbered by position in the figure
TMAX = 11.3
TMIN = -1.9      # room on the left for the fiber with the gratings
FIB_T = -1.0     # where the fiber is drawn


def tint(col, w):
    """Blend a colour with light grey: w = 1 keeps the colour, w = 0 is grey."""
    c = np.array(to_rgb(col)); g = np.array([0.78, 0.78, 0.78])
    return tuple(w * c + (1 - w) * g)


# ghosts take the colour of the grating of their first reflection, greyed
GH = {g: tint(COL[g], 0.45) for g in COL}    # lines, bars, dots
GHT = {g: tint(COL[g], 0.7) for g in COL}    # text


def z(g):
    return T[g] / 2.0


def path_points(seq):
    """(t, z) vertices of the path that reflects at the gratings in seq,
    starting from the launch at (0, 0) and ending at the photodiode."""
    pts = [(0.0, 0.0)]
    t, zz = 0.0, 0.0
    for g in seq:
        t += abs(z(g) - zz)
        zz = z(g)
        pts.append((t, zz))
    t += zz
    pts.append((t, 0.0))
    return np.array(pts)


GHOST_SEQS = [('c', 'b', 'c'), ('a', 'b', 'c'), ('c', 'b', 'a'), ('a', 'c', 'a'), ('a', 'b', 'a')]
GHOST_LABELS = {2 * TC - TB: r'$2\tau_2{-}\tau_1$',
                TA - TB + TC: r'$\tau_3{-}\tau_1{+}\tau_2$',
                2 * TA - TC: r'$2\tau_3{-}\tau_2$',
                2 * TA - TB: r'$2\tau_3{-}\tau_1$'}

fig = plt.figure(figsize=(3.45, 5.0))
ax = fig.add_axes([0.105, 0.665, 0.88, 0.295])
cx = fig.add_axes([0.105, 0.365, 0.88, 0.205])
arrivals = fig.add_axes([0.105, 0.582, 0.88, 0.05])
dx = fig.add_axes([0.105, 0.08, 0.88, 0.205])
DIRECT, GHOST = FS.BLUE, FS.VERM

LW_IN, LW_RET, LW_GH, LW_GHF = 2.2, 1.5, 1.4, 0.7


def seg_arrow(a_, p0, p1, color, lw, ms=7, frac=0.55):
    """A line segment with an arrowhead at fraction frac of its length."""
    a_.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=lw, solid_capstyle='round', zorder=3)
    q = (p0[0] + frac * (p1[0] - p0[0]), p0[1] + frac * (p1[1] - p0[1]))
    a_.add_patch(FancyArrowPatch(p0, q, arrowstyle='-|>', mutation_scale=ms, color=color,
                                 lw=0, shrinkA=0, shrinkB=0, zorder=4))


# ---------------- (a) space-time diagram -------------------------------------
# Draw the direct returns and one complete ghost path. The remaining paths
# are counted in the separate arrival strip above panel (b).
ZMAX, diagram_end = z('a') + 0.40, 8.5
ax.plot([FIB_T, FIB_T], [0, ZMAX - 0.1], color='0.65', lw=1.2, zorder=1)
ax.add_patch(Rectangle((FIB_T - 0.17, -0.08), 0.34, 0.16,
                       facecolor='0.3', edgecolor='none'))
ax.text(FIB_T - 0.28, 0, 'PD', ha='right', va='center', fontsize=6.2,
         color='0.30')
ax.text(FIB_T - 0.1, ZMAX + 0.02, 'Position $z$', ha='center',
         va='bottom', fontsize=6.0, color='0.35')
for g in 'bca':
    zg = z(g)
    ax.plot([FIB_T + 0.2, diagram_end], [zg, zg], color='0.82',
             lw=0.5, ls=(0, (2, 2)), zorder=0)
    ax.add_patch(Rectangle((FIB_T - 0.18, zg - 0.09), 0.36, 0.18,
                           facecolor='white', edgecolor=DIRECT, lw=0.7, zorder=4))
    for xg in np.linspace(FIB_T - 0.12, FIB_T + 0.12, 4):
        ax.plot([xg, xg], [zg - 0.065, zg + 0.065],
                 color=DIRECT, lw=0.55, zorder=5)
    ax.text(FIB_T - 0.30, zg, 'FBG ' + NAME[g], ha='right', va='center',
             fontsize=6.2, color='0.25')

seg_arrow(ax, (0, 0), (z('a'), z('a')), '0.40', 1.1, ms=5.5)
ax.text(0.91, 1.56, 'Launch', color='0.4', fontsize=6.0,
         rotation=0, ha='center')
for g in 'bca':
    direct_path = path_points((g,))
    seg_arrow(ax, tuple(direct_path[1]), tuple(direct_path[2]),
              DIRECT, 0.9, ms=5, frac=0.70)
    ax.plot(T[g], 0, 's', color=DIRECT, ms=3.2, zorder=6)
    ax.text(T[g], -0.15, r'$\tau_%s$' % NAME[g], color=DIRECT,
             ha='center', va='top', fontsize=6.3)

ghost_path = path_points(('a', 'b', 'c'))
for i in range(1, len(ghost_path) - 1):
    seg_arrow(ax, tuple(ghost_path[i]), tuple(ghost_path[i + 1]),
              GHOST, 1.45, ms=5.6, frac=0.65)
ax.plot(ghost_path[1:-1, 0], ghost_path[1:-1, 1], 'o', color=GHOST,
         ms=3.0, mec='white', mew=0.5, zorder=6)
ax.plot(ghost_path[-1, 0], 0, 'o', color=GHOST, ms=3.5, zorder=6)
ax.text(ghost_path[-1, 0], -0.15, r'$\tau_g$', ha='center', va='top',
         fontsize=6.3, color=GHOST)
ax.text(8.32, 2.66, r'Ghost $3\to1\to2$', ha='right',
         fontsize=6.6, color=GHOST)
ax.text(8.32, 2.24, r'$\tau_g=\tau_3-\tau_1+\tau_2$', ha='right',
         fontsize=6.2, color=GHOST)
ax.text(8.32, 1.92, r'$P_g\propto R_3R_1R_2$', ha='right',
         fontsize=6.2, color=GHOST)
ax.text(1.55, 0.25, 'Direct', fontsize=5.8, color=DIRECT, ha='center')
ax.add_patch(FancyArrowPatch((0, 0), (diagram_end, 0), arrowstyle='-|>',
                             mutation_scale=6, lw=0.7, color='0.3', zorder=1))
ax.text(diagram_end, -0.15, r'Delay $\tau$', ha='right', va='top',
         fontsize=6.2, color='0.25')
ax.set_xlim(TMIN, diagram_end)
ax.set_ylim(-0.40, ZMAX)
ax.axis('off')

# ---------------- (b) delay profile with the arrivals as bars ---------------
# Receiver chain as in Section IV: launched power P0 into a link of loss ALPHA,
# every path up to the third order, photodiode of responsivity RESP with the
# NEP of the link budget, shot noise of the photocurrent, laser RIN, a 4th-order
# Bessel low-pass at 0.75 B, a 12-bit ADC at NS samples per chip, record mean
# removed, periodic correlation with the bipolar replica, constant offset of
# the periodic correlation (-1/N per return) subtracted.
from scipy.signal import bessel, filtfilt
C_LIGHT, NG = 2.998e8, 1.468
CHIP_RATE, NBITS, R = 25e6, 7, 0.10
SPC_FINE, NS = 32, 8                                 # fine grid and ADC samples per chip
Z = {'b': 8.0, 'c': 19.2, 'a': 40.0}                 # metres, the ratios of (a)
M_PER_CHIP = C_LIGHT / NG / CHIP_RATE / 2.0
P0, ALPHA, NEP, RESP = 1e-3, 10 ** (-4.0 / 10), 0.5e-12, 0.9   # W, link loss, W/sqrt(Hz), A/W
RIN_DB, ADC_BITS = -130.0, 12
Q_E = 1.602e-19
RNG = np.random.default_rng(3)


def paths(Zd):
    """[(sequence of gratings, amplitude, equivalent position), ...] up to the third order."""
    order = sorted(Zd, key=Zd.get)

    def trans(z0, z1):
        t = 1.0
        for g in order:
            if min(z0, z1) < Zd[g] < max(z0, z1):
                t *= (1.0 - R)
        return t

    def one(seq):
        amp, zz, length = 1.0, 0.0, 0.0
        for g in seq:
            amp *= trans(zz, Zd[g]) * R
            length += abs(Zd[g] - zz)
            zz = Zd[g]
        amp *= trans(zz, 0.0)
        return amp, (length + zz) / 2.0

    out = [((g,),) + one((g,)) for g in order]
    for g1 in order:
        for g2 in order:
            for g3 in order:
                if Zd[g2] < Zd[g1] and Zd[g2] < Zd[g3]:
                    out.append(((g1, g2, g3),) + one((g1, g2, g3)))
    return out


def delay_profile(Zd, noise=True):
    code01 = C._mls01(NBITS)
    n = code01.size
    fs = CHIP_RATE * SPC_FINE
    tx = np.repeat(code01.astype(float), SPC_FINE)
    # optical power at the photodiode, one code period, every path
    popt = np.zeros_like(tx)
    for seq, amp, zpos in paths(Zd):
        popt += amp * np.roll(tx, int(round(zpos / M_PER_CHIP * SPC_FINE)))
    popt *= P0 * ALPHA
    if noise:
        popt *= 1.0 + RNG.normal(0.0, np.sqrt(10 ** (RIN_DB / 10) * fs / 2), popt.size)   # laser RIN
    i_pd = RESP * popt
    if noise:
        i_pd += RNG.normal(0.0, np.sqrt(2 * Q_E * RESP * popt.mean() * fs / 2), i_pd.size)   # shot noise
        i_pd += RNG.normal(0.0, RESP * NEP * np.sqrt(fs / 2), i_pd.size)                     # receiver NEP
    b_, a_ = bessel(4, 0.75 * CHIP_RATE / (fs / 2), norm='mag')
    i_pd = filtfilt(b_, a_, i_pd)                    # receiver low-pass, zero phase
    dec = SPC_FINE // NS
    rec = i_pd[dec // 2::dec]                        # ADC samples, NS per chip
    fsr = 1.2 * rec.max()
    rec = np.round(rec / fsr * 2 ** (ADC_BITS - 1)) / 2 ** (ADC_BITS - 1) * fsr   # quantization
    replica = np.repeat(C._to_pm1(code01).astype(float), NS)
    rec = rec - rec.mean()
    corr = np.fft.ifft(np.fft.fft(rec) * np.conj(np.fft.fft(replica))).real
    corr /= -(n * NS / 2.0) * RESP * P0 * ALPHA     # a unit return gives a peak of +1
    corr -= np.median(corr)                          # constant offset of the periodic correlation
    return np.arange(rec.size) / NS * M_PER_CHIP, corr, n


zaxis, corr, n = delay_profile(Z)
P = paths(Z)
ref = P[0][1]                                      # ideal amplitude of the first direct return
ref_meas = corr[(zaxis > 6.0) & (zaxis < 10.0)].max()   # its peak after the receiver chain
print('first direct return after the receiver chain: %.3f of ideal' % (ref_meas / ref))
XMAX = 80.0
m = zaxis <= XMAX
# The arrival strip encodes path count and delay, never optical amplitude.
arrivals.set_xlim(0, XMAX)
arrivals.set_ylim(0, 2.7)
arrivals.axis('off')
arrivals.plot([0, XMAX], [0.65, 0.65], color='0.85', lw=0.5)
ghost_count = {}
for seq, amp, zpos in P:
    if len(seq) == 1:
        arrivals.plot(zpos, 0.65, 's', color=DIRECT, ms=3.4)
        arrivals.text(zpos, 1.40, r'$\tau_%s$' % NAME[seq[0]], color=DIRECT,
                       ha='center', va='bottom', fontsize=6.1)
    else:
        stack = ghost_count.get(zpos, 0)
        arrivals.plot(zpos, 0.65 + 0.72 * stack, 'o', color=GHOST, ms=3.1)
        ghost_count[zpos] = stack + 1
arrivals.text(0, 2.50, 'Arrival times', fontsize=5.8, color='0.40',
               ha='left', va='bottom')
cx.plot(zaxis[m], corr[m] / ref_meas, color='0.22', lw=0.95, zorder=3)
cx.set_xlim(0, XMAX)
cx.set_ylim(0, 1.12)
cx.set_xticks([0, 20, 40, 60, 80])
cx.set_yticks([0, 0.5, 1.0])
cx.tick_params(labelsize=6.2, length=2.4)
cx.set_xlabel('Equivalent position (m)', fontsize=6.8, labelpad=1.5)
cx.set_ylabel('Normalized correlation', fontsize=6.8, labelpad=1.5)

zoom = cx.inset_axes([0.62, 0.44, 0.35, 0.43])
late = (zaxis >= 46) & (zaxis <= 77)
zoom.plot(zaxis[late], corr[late] / ref_meas, color='0.22', lw=0.8)
for zpos in sorted(ghost_count):
    if zpos > 46:
        zoom.axvline(zpos, color=GHOST, ls=(0, (2, 2)), lw=0.55, zorder=0)
zoom.set_xlim(46, 77)
zoom.set_ylim(0, 0.021)
zoom.set_xticks([50, 60, 70])
zoom.set_yticks([0, 0.01, 0.02])
zoom.tick_params(labelsize=4.8, length=1.7, pad=1.2)
zoom.minorticks_off()
for spine in zoom.spines.values():
    spine.set_linewidth(0.65)
zoom.set_title('Late echoes', fontsize=5.8, pad=2.2, loc='center')

# ---------------- (c) what deshadowing does to a ghost -----------------------
# Spectral model of s41: eight gratings at R = 10 %, uniform 4-m spacing,
# lines spread over the band, every third-order ghost summed under each
# grating, the -1/N leakage of the other lines, then the sequential
# correction (14). Shown for the last grating, which carries the most ghosts.
K8, R8, FWHM8, N8, B8 = 8, 0.10, 250.0, 127, 100e6
NU0 = np.linspace(-175.0, 175.0, K8)


def sline(nu, nu0):
    return np.exp(-0.5 * ((nu - nu0) / (FWHM8 / 2.35482)) ** 2)


def spectral(zpos, ghosts=True):
    nu = np.linspace(-520.0, 520.0, 1041)
    tau = 2.0 * NG * zpos / C_LIGHT * B8
    shapes = np.array([sline(nu, c) for c in NU0])
    trans = np.ones((K8, nu.size))
    for k in range(1, K8):
        trans[k] = trans[k - 1] * (1.0 - R8 * shapes[k - 1]) ** 2
    A = R8 * shapes * trans
    ghost = np.zeros_like(A)
    if ghosts:
        for b in range(K8):
            for a in range(K8):
                for c in range(K8):
                    if b < a and b < c:
                        tg = tau[a] - tau[b] + tau[c]
                        for k in range(K8):
                            w = 1.0 - abs(tau[k] - tg)
                            if w > 0.02:
                                ghost[k] += w * R8 ** 3 * shapes[a] * shapes[b] * shapes[c]
    meas = np.empty_like(A)
    for k in range(K8):
        meas[k] = A[k] + ghost[k] - (1.0 / N8) * (A.sum(axis=0) - A[k])
    est = np.ones_like(nu)
    corr = np.empty_like(A)
    for k in range(K8):
        corr[k] = meas[k] / np.maximum(est, 0.05)
        est = est * (1.0 - np.clip(corr[k], 0.0, 0.99)) ** 2
    return nu, shapes, meas, corr


z8 = 4.0 * np.arange(1, K8 + 1)
nu8, shapes8, meas8, corr8 = spectral(z8)
_, _, _, corr8_0 = spectral(z8, ghosts=False)
k8 = K8 - 1
dx.plot(nu8, shapes8[k8], color='0.6', lw=1.8, label='line of grating 8', zorder=2)
dx.plot(nu8, meas8[k8] / R8, color=COL['b'], lw=1.0, label='measured', zorder=3)
dx.plot(nu8, corr8[k8] / R8, color=FS.BLUE, lw=1.0, label='deshadowed', zorder=4)
dx.plot(nu8, (corr8[k8] - corr8_0[k8]) / R8, color=COL['a'], lw=0.9, ls='--', label='ghosts left after deshadowing', zorder=3)
dx.set_xlim(-500, 500)
dx.set_ylim(0, 1.12)
dx.set_yticks([0, 0.5, 1.0])
dx.tick_params(labelsize=6.2, length=2.4)
dx.set_xlabel('wavelength offset (pm)', fontsize=6.8, labelpad=1.5)
dx.set_ylabel('reflectance / $R$', fontsize=6.8, labelpad=1.5)
dx.legend(fontsize=5.8, loc='upper left', frameon=True, handlelength=1.8, borderpad=0.4, labelspacing=0.25)
dx.text(0.985, 0.95, '8 gratings, uniform 4 m, $R=10\%$', transform=dx.transAxes, ha='right', va='top', fontsize=5.8, color=GREY)

for letter, ypos in (('a', 0.976), ('b', 0.644), ('c', 0.304)):
    fig.text(0.028, ypos, letter, fontsize=9, fontweight='bold', va='bottom')

os.makedirs('figs', exist_ok=True)
fig.savefig('figs/fig_s50_ghostpath.pdf')
fig.savefig('figs/fig_s50_ghostpath.png', dpi=300)
print('saved figs/fig_s50_ghostpath.pdf')
for seq, amp, zpos in P:
    print('  %-10s amp %.2e  z %5.1f m' % ('(' + ','.join(seq) + ')', amp, zpos))
