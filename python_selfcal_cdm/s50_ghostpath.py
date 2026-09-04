"""s50_ghostpath.py - every arrival at the photodiode for three gratings.

Replaces the TikZ drawing fig_ghost.tex with a matplotlib figure in the house
style (figstyle), so that fonts, colours and panel letters match the other
figures of the paper.

  (a) Space-time diagram. Horizontal axis: time since the launch, in units of
      the delay. Vertical axis: position along the fibre, with the fibre
      itself, its three gratings 1, 2, 3 and the photodiode drawn on the
      left. The launched code climbs the diagonal; every reflection descends
      back to z = 0 and lands on the time axis at its delay. The three direct
      returns and the third-order path (3,1,2) are drawn in full, the other
      four third-order paths as thin lines. Every ghost carries the colour of
      the grating of its first reflection, greyed.
  (b) The delay profile of the same three gratings in the time-domain model:
      z = 8, 19.2, 40 m (the ratios of (a)), R = 10 %, unipolar m-sequence
      N = 127 at 25 Mchip/s (Table III), every path up to the third order,
      record mean removed, correlation with the bipolar replica, constant
      offset of the periodic correlation subtracted, linear scale. Behind the
      curve, a coloured bar marks every arrival with its delay: three direct
      returns of power ~ R and five third-order paths of power ~ R^3, bar
      heights not to scale. The spacing is non-uniform on purpose so that the
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

fig = plt.figure(figsize=(3.45, 3.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1.5, 1.0],
                      hspace=0.2, left=0.075, right=0.985, bottom=0.11, top=0.95)
ax = fig.add_subplot(gs[0, 0])
cx = fig.add_subplot(gs[1, 0])

LW_IN, LW_RET, LW_GH, LW_GHF = 2.2, 1.5, 1.4, 0.7


def seg_arrow(a_, p0, p1, color, lw, ms=7, frac=0.55):
    """A line segment with an arrowhead at fraction frac of its length."""
    a_.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=lw, solid_capstyle='round', zorder=3)
    q = (p0[0] + frac * (p1[0] - p0[0]), p0[1] + frac * (p1[1] - p0[1]))
    a_.add_patch(FancyArrowPatch(p0, q, arrowstyle='-|>', mutation_scale=ms, color=color,
                                 lw=0, shrinkA=0, shrinkB=0, zorder=4))


# ---------------- (a) space-time diagram -------------------------------------
ZMAX = z('a') + 0.55
# the fiber itself, drawn vertically on the left, with the photodiode at the
# bottom and a grating symbol (short stripes) at every z_k
ax.plot([FIB_T, FIB_T], [0.0, ZMAX - 0.12], color='0.72', lw=3.4, solid_capstyle='butt', zorder=2)
ax.plot([FIB_T, FIB_T], [0.0, ZMAX - 0.12], color='0.5', lw=0.7, solid_capstyle='butt', zorder=2)
ax.add_patch(FancyArrowPatch((FIB_T, ZMAX - 0.4), (FIB_T, ZMAX - 0.02), arrowstyle='-|>', mutation_scale=7,
                             color='0.5', lw=0.7, shrinkA=0, shrinkB=0, zorder=3))
ax.text(FIB_T - 0.32, ZMAX - 0.05, 'position $z$', color=GREY, ha='right', va='top', fontsize=6.2)
ax.add_patch(Rectangle((FIB_T - 0.25, -0.1), 0.5, 0.2, facecolor='0.35', edgecolor='none', zorder=4))
ax.text(FIB_T - 0.42, 0.0, 'PD', color=GREY, ha='right', va='center', fontsize=6.6)
for g in 'bca':
    zg = z(g)
    ax.add_patch(Rectangle((FIB_T - 0.24, zg - 0.1), 0.48, 0.2, facecolor='white', edgecolor=COL[g], lw=0.7, zorder=3))
    for k in range(5):
        xk = FIB_T - 0.16 + 0.08 * k
        ax.plot([xk, xk], [zg - 0.07, zg + 0.07], color=COL[g], lw=0.7, zorder=4)
    ax.plot([FIB_T + 0.26, TMAX], [zg, zg], color=COL[g], lw=0.9, ls=(0, (1.2, 1.8)), zorder=1)
    ax.text(FIB_T - 0.32, zg, 'FBG %s' % NAME[g], color=COL[g], ha='right', va='center', fontsize=6.8)

# the launched code climbs the diagonal past every grating
seg_arrow(ax, (0, 0), (z('a'), z('a')), '0.5', LW_IN, ms=8, frac=0.45)
ax.plot([z('a'), z('a') + 0.45], [z('a'), z('a') + 0.45], color='0.5', lw=LW_IN, ls=(0, (1.2, 1.4)), zorder=2)
ax.text(1.85, 2.12, 'launched code', color=GREY, ha='right', va='center', fontsize=6.2)

# direct returns
for g in 'bca':
    p = path_points((g,))
    seg_arrow(ax, tuple(p[1]), tuple(p[2]), COL[g], LW_RET, frac=0.6)
    ax.text(p[1][0] + 0.12, p[1][1] - 0.02, r'$\times R_%s$' % NAME[g], color=COL[g], ha='left', va='top', fontsize=6.6)

# the other third-order paths, thin
for seq in GHOST_SEQS:
    if seq == ('a', 'b', 'c'):
        continue
    p = path_points(seq)
    ax.plot(p[1:, 0], p[1:, 1], color=GH[seq[0]], lw=LW_GHF + 0.2, zorder=2)

# the path (3,1,2) in full
p = path_points(('a', 'b', 'c'))
seg_arrow(ax, tuple(p[1]), tuple(p[2]), GH['a'], LW_GH, ms=6, frac=0.55)
seg_arrow(ax, tuple(p[2]), tuple(p[3]), GH['a'], LW_GH, ms=6, frac=0.55)
seg_arrow(ax, tuple(p[3]), tuple(p[4]), GH['a'], LW_GH, ms=6, frac=0.55)
ax.text(p[2][0] + 0.1, p[2][1] + 0.05, r'$\times R_1$', color=COL['b'], ha='left', va='bottom', fontsize=6.6)
ax.text(p[3][0] + 0.12, p[3][1] - 0.02, r'$\times R_2$', color=COL['c'], ha='left', va='top', fontsize=6.6)
ax.text(TMAX - 0.1, 2.05, r'path $(3,1,2)$: $\tau_3{-}\tau_1{+}\tau_2$', color=GHT['a'], ha='right', va='bottom', fontsize=6.4)
ax.text(TMAX - 0.1, 1.82, r'power $\propto R_3R_1R_2$', color=GHT['a'], ha='right', va='bottom', fontsize=6.4)

# arrivals at the photodiode: dots on the time axis
for g in 'bca':
    ax.plot(T[g], 0, 'o', color=COL[g], ms=3.6, mec='white', mew=0.5, zorder=5)
for seq in GHOST_SEQS:
    ax.plot(path_points(seq)[-1, 0], 0, 'o', color=GH[seq[0]], ms=3.0, mec='white', mew=0.5, zorder=5)
ax.text(TMAX - 0.7, 0.12, r'time $\tau$', color='0.25', ha='right', va='bottom', fontsize=7)
ax.add_patch(FancyArrowPatch((TMAX - 0.6, 0), (TMAX, 0), arrowstyle='-|>', mutation_scale=7, color='0.3',
                             lw=0.9, shrinkA=0, shrinkB=0, zorder=4, clip_on=False))

ax.set_xlim(TMIN, TMAX)
ax.set_ylim(-0.14, ZMAX)
ax.set_yticks([]); ax.set_xticks([])
for s in ('top', 'right', 'left'):
    ax.spines[s].set_visible(False)
ax.spines['bottom'].set_color('0.3')
ax.spines['bottom'].set_position(('data', 0))
ax.spines['bottom'].set_bounds(0, TMAX)

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
# every arrival as a bar, heights schematic (direct returns 1, ghosts 0.33)
H_DIR, H_GH = 1.0, 0.33
LABEL = {('c', 'b', 'c'): (r'$2\tau_2{-}\tau_1$', 0), ('a', 'b', 'c'): (r'$\tau_3{-}\tau_1{+}\tau_2$', 1),
         ('a', 'c', 'a'): (r'$2\tau_3{-}\tau_2$', 0), ('a', 'b', 'a'): (r'$2\tau_3{-}\tau_1$', 1)}
for seq, amp, zpos in P:
    if len(seq) == 1:
        cx.plot([zpos, zpos], [0, H_DIR], color=COL[seq[0]], lw=2.4, solid_capstyle='butt', zorder=2)
        cx.text(zpos, H_DIR + 0.05, r'$\tau_%s$' % NAME[seq[0]], color=COL[seq[0]], ha='center', va='bottom', fontsize=6.4)
    else:
        off = 0.45 if seq == ('a', 'b', 'c') else (-0.45 if seq == ('c', 'b', 'a') else 0.0)
        cx.plot([zpos + off, zpos + off], [0, H_GH], color=GH[seq[0]], lw=1.8, solid_capstyle='butt', zorder=2)
        if seq in LABEL:
            lab, row = LABEL[seq]
            col = '0.4' if seq == ('a', 'b', 'c') else GHT[seq[0]]
            cx.text(zpos, H_GH + 0.06 + 0.17 * row, lab, color=col, ha='center', va='bottom', fontsize=5.8)
cx.plot(zaxis[m], corr[m] / ref_meas, color='0.25', lw=0.8, zorder=3)
cx.set_xlim(0, XMAX)
cx.set_ylim(0, 1.2)
cx.set_xticks([0, 20, 40, 60, 80])
cx.set_yticks([0, 0.5, 1.0])
cx.tick_params(labelsize=6.2, length=2.4)
cx.set_xlabel('position (m)', fontsize=6.8, labelpad=1.5)
cx.set_ylabel('correlation, norm.', fontsize=6.8, labelpad=1.5)
cx.text(0.985, 0.95, r'$R=10\%$, $N=127$', transform=cx.transAxes, ha='right', va='top', fontsize=6.0)
cx.text(0.985, 0.83, 'bar heights not to scale', transform=cx.transAxes, ha='right', va='top', fontsize=5.8, color=GREY)
for a_, let in ((ax, 'a'), (cx, 'b')):
    a_.text(-0.02, 1.02, let, transform=a_.transAxes, fontsize=9, fontweight='bold',
            ha='right', va='bottom')

os.makedirs('figs', exist_ok=True)
fig.savefig('figs/fig_s50_ghostpath.pdf')
fig.savefig('figs/fig_s50_ghostpath.png', dpi=220)
print('saved figs/fig_s50_ghostpath.pdf')
for seq, amp, zpos in P:
    print('  %-10s amp %.2e  z %5.1f m' % ('(' + ','.join(seq) + ')', amp, zpos))
