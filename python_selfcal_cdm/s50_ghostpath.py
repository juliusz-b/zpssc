"""s50_ghostpath.py - every arrival at the photodiode for three gratings.

Replaces the TikZ drawing fig_ghost.tex with a matplotlib figure in the house
style (figstyle), so that fonts, colours and panel letters match the other
figures of the paper.

  (a) Space-time diagram. Horizontal axis: time since the launch, in units of
      the delay. Vertical axis: position along the fibre, with the photodiode
      at z = 0 and the three gratings a, b, c at z_k = c tau_k / (2 n_g). The
      launched code climbs the diagonal; every reflection descends back to
      z = 0 and lands on the time axis at its delay. The three direct returns
      and the third-order path (a,b,c) are drawn in full, the other four
      third-order paths as thin lines.
  (b) The same time axis with every arrival: three direct returns of power
      ~ R and five third-order paths of power ~ R^3, labelled with delay and
      reflectivity product. The spacing is non-uniform on purpose so that the
      ghosts do not collide, except the pair (a,b,c)/(c,b,a), which share one
      delay for any spacing.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import figstyle as FS

FS.apply()

# ---------------------------------------------------------------------------
# geometry: time in delay units, position z = tau / 2 (c / n_g = 1)
# ---------------------------------------------------------------------------
TB, TC, TA = 1.0, 2.4, 5.0
T = {'a': TA, 'b': TB, 'c': TC}
COL = {'a': FS.VERM, 'b': FS.ORANGE, 'c': FS.GREEN}
PUR = FS.PURPLE
GREY = '0.45'
TMAX = 11.3


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
GHOST_LABELS = {2 * TC - TB: (r'$2\tau_c{-}\tau_b$', r'$R_bR_c^{2}$'),
                TA - TB + TC: (r'$\tau_a{-}\tau_b{+}\tau_c$', r'$2R_aR_bR_c$'),
                2 * TA - TC: (r'$2\tau_a{-}\tau_c$', r'$R_a^{2}R_c$'),
                2 * TA - TB: (r'$2\tau_a{-}\tau_b$', r'$R_a^{2}R_b$')}

fig = plt.figure(figsize=(3.45, 3.4))
gs = fig.add_gridspec(2, 1, height_ratios=[1.55, 1.0], hspace=0.06,
                      left=0.115, right=0.985, bottom=0.05, top=0.94)
ax = fig.add_subplot(gs[0, 0])
bx = fig.add_subplot(gs[1, 0], sharex=ax)

LW_IN, LW_RET, LW_GH, LW_GHF = 2.2, 1.5, 1.4, 0.7


def seg_arrow(a_, p0, p1, color, lw, ms=7, frac=0.55):
    """A line segment with an arrowhead at fraction frac of its length."""
    a_.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=lw, solid_capstyle='round', zorder=3)
    q = (p0[0] + frac * (p1[0] - p0[0]), p0[1] + frac * (p1[1] - p0[1]))
    a_.add_patch(FancyArrowPatch(p0, q, arrowstyle='-|>', mutation_scale=ms, color=color,
                                 lw=0, shrinkA=0, shrinkB=0, zorder=4))


# ---------------- (a) space-time diagram -------------------------------------
ZMAX = z('a') + 0.55
# gratings: horizontal dotted lines with labels on the left
for g in 'bca':
    ax.plot([0, TMAX], [z(g), z(g)], color=COL[g], lw=0.9, ls=(0, (1.2, 1.8)), zorder=1)
    ax.text(-0.12, z(g), r'$%s$, $R_%s$' % (g, g), color=COL[g], ha='right', va='center', fontsize=7)
ax.text(-0.12, ZMAX - 0.05, 'position $z$', color=GREY, ha='right', va='top', fontsize=6.2)
ax.text(-0.12, 0.0, 'PD', color=GREY, ha='right', va='center', fontsize=6.6)

# the launched code climbs the diagonal past every grating
seg_arrow(ax, (0, 0), (z('a'), z('a')), '0.5', LW_IN, ms=8, frac=0.45)
ax.plot([z('a'), z('a') + 0.45], [z('a'), z('a') + 0.45], color='0.5', lw=LW_IN, ls=(0, (1.2, 1.4)), zorder=2)
ax.text(1.8, 2.12, 'launched code', color=GREY, ha='right', va='center', fontsize=6.2)

# direct returns
for g in 'bca':
    p = path_points((g,))
    seg_arrow(ax, tuple(p[1]), tuple(p[2]), COL[g], LW_RET, frac=0.6)
    ax.text(p[1][0] + 0.12, p[1][1] - 0.02, r'$\times R_%s$' % g, color=COL[g], ha='left', va='top', fontsize=6.6)

# the other third-order paths, thin
for seq in GHOST_SEQS:
    if seq == ('a', 'b', 'c'):
        continue
    p = path_points(seq)
    ax.plot(p[1:, 0], p[1:, 1], color=PUR, lw=LW_GHF, alpha=0.5, zorder=2)

# the path (a,b,c) in full
p = path_points(('a', 'b', 'c'))
seg_arrow(ax, tuple(p[1]), tuple(p[2]), PUR, LW_GH, ms=6, frac=0.55)
seg_arrow(ax, tuple(p[2]), tuple(p[3]), PUR, LW_GH, ms=6, frac=0.55)
seg_arrow(ax, tuple(p[3]), tuple(p[4]), PUR, LW_GH, ms=6, frac=0.55)
ax.text(p[2][0] + 0.1, p[2][1] + 0.05, r'$\times R_b$', color=COL['b'], ha='left', va='bottom', fontsize=6.6)
ax.text(p[3][0] + 0.12, p[3][1] - 0.02, r'$\times R_c$', color=COL['c'], ha='left', va='top', fontsize=6.6)
ax.text(TMAX - 0.1, 2.05, r'path $(a,b,c)$: $\tau_a{-}\tau_b{+}\tau_c$', color=PUR, ha='right', va='bottom', fontsize=6.4)
ax.text(TMAX - 0.1, 1.82, r'power $\propto R_aR_bR_c$', color=PUR, ha='right', va='bottom', fontsize=6.4)

# arrivals at the photodiode: dots on the time axis
for g in 'bca':
    ax.plot(T[g], 0, 'o', color=COL[g], ms=3.6, mec='white', mew=0.5, zorder=5)
for seq in GHOST_SEQS:
    ax.plot(path_points(seq)[-1, 0], 0, 'o', color=PUR, ms=3.0, mec='white', mew=0.5, zorder=5)

ax.set_xlim(-0.05, TMAX)
ax.set_ylim(-0.12, ZMAX)
ax.set_yticks([]); ax.set_xticks([])
for s in ('top', 'right', 'left'):
    ax.spines[s].set_visible(False)
ax.spines['bottom'].set_color('0.3')
ax.spines['bottom'].set_position(('data', 0))

# ---------------- (b) arrivals on the same time axis ------------------------
H_DIR, H_GH = 1.0, 0.33
for g in 'bca':
    bx.plot([T[g], T[g]], [0, H_DIR], color=COL[g], lw=3.0, solid_capstyle='butt', zorder=3)
    bx.text(T[g], H_DIR + 0.06, r'$R_%s$' % g, color=COL[g], ha='center', va='bottom', fontsize=7)
    bx.text(T[g], -0.1, r'$\tau_%s$' % g, color=COL[g], ha='center', va='top', fontsize=7)
for tg, (lab_t, lab_r) in GHOST_LABELS.items():
    double = abs(tg - (TA - TB + TC)) < 1e-9
    xs = (tg - 0.09, tg + 0.09) if double else (tg,)
    for x_ in xs:
        bx.plot([x_, x_], [0, H_GH], color=PUR, lw=2.2, solid_capstyle='butt', zorder=3)
    bx.text(tg, H_GH + 0.06 + (0.22 if double else 0), lab_r, color=PUR, ha='center', va='bottom', fontsize=6.6)
    bx.text(tg, -0.1 - (0.24 if double else 0), lab_t, color=PUR, ha='center', va='top', fontsize=6.4)
# leaders from the arrival dots to the bars
for g in 'bca':
    bx.plot([T[g], T[g]], [H_DIR + 0.02, 1.62], color=COL[g], lw=0.7, ls=(0, (1, 1.6)), alpha=0.8, zorder=0)
for tg in GHOST_LABELS:
    bx.plot([tg, tg], [H_GH + (0.36 if abs(tg - (TA - TB + TC)) < 1e-9 else 0.3), 1.62],
            color=PUR, lw=0.7, ls=(0, (1, 1.6)), alpha=0.6, zorder=0)
bx.text(TMAX - 0.05, 0.12, r'time $\tau$', color='0.25', ha='right', va='bottom', fontsize=7)
for a_ in (ax, bx):
    a_.add_patch(FancyArrowPatch((TMAX - 0.6, 0), (TMAX, 0), arrowstyle='-|>', mutation_scale=7, color='0.3', lw=0.9, shrinkA=0, shrinkB=0, zorder=4, clip_on=False))
bx.text(0.0, 1.55, 'arrivals at the PD, heights not to scale', color=GREY, ha='left', va='top', fontsize=6.2)
bx.set_ylim(-0.6, 1.62)
bx.set_yticks([]); bx.set_xticks([])
for s in ('top', 'right', 'left'):
    bx.spines[s].set_visible(False)
bx.spines['bottom'].set_color('0.3')
bx.spines['bottom'].set_position(('data', 0))

for a_, let in ((ax, 'a'), (bx, 'b')):
    a_.text(-0.02, 1.02, let, transform=a_.transAxes, fontsize=9, fontweight='bold', ha='right', va='bottom')

os.makedirs('figs', exist_ok=True)
fig.savefig('figs/fig_s50_ghostpath.pdf')
fig.savefig('figs/fig_s50_ghostpath.png', dpi=220)
print('saved figs/fig_s50_ghostpath.pdf')
