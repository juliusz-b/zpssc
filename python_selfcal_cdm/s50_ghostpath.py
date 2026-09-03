"""s50_ghostpath.py - every arrival at the photodiode for three gratings.

Replaces the TikZ drawing fig_ghost.tex with a matplotlib figure in the house
style (figstyle), so that fonts, colours and panel letters match the other
figures of the paper.

  (a) The fibre with three gratings a, b, c drawn at x proportional to their
      delay, the launched code, the three direct returns on their own lanes,
      and one of the five third-order paths, (a,b,c): it branches off the
      return of a at b and is re-reflected at c.
  (b) The delay axis, aligned with the fibre above, with every arrival at the
      photodiode: three direct returns of power ~ R and five third-order paths
      of power ~ R^3, each labelled with its delay and its reflectivity product.
      The spacing is non-uniform on purpose so that the ghosts do not collide,
      except the pair (a,b,c)/(c,b,a), which share a delay for any spacing.
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
# geometry, x in delay units
# ---------------------------------------------------------------------------
TB, TC, TA = 1.0, 2.4, 5.0            # delays of b, c, a
X0 = -0.55                            # photodiode
XEND = 6.0                            # end of the drawn fibre
GHOSTS = [(2 * TC - TB, r'$2\tau_c{-}\tau_b$', r'$R_bR_c^{2}$', 1),
          (TA - TB + TC, r'$\tau_a{-}\tau_b{+}\tau_c$', r'$2R_aR_bR_c$', 2),
          (2 * TA - TC, r'$2\tau_a{-}\tau_c$', r'$R_a^{2}R_c$', 1),
          (2 * TA - TB, r'$2\tau_a{-}\tau_b$', r'$R_a^{2}R_b$', 1)]
GR = [(TB, 'b', FS.ORANGE), (TC, 'c', FS.GREEN), (TA, 'a', FS.VERM)]
PUR = FS.PURPLE
GREY = '0.45'

# lanes (y) of panel (a)
Y_FIB = 0.0
Y_IN = -0.55
Y_RB, Y_RC, Y_RA = -0.95, -1.35, -1.75      # direct returns of b, c, a
Y_G1, Y_G2 = -2.15, -2.55                    # ghost: forward leg, return leg
R_HOOK = 0.18                                # hook radius in x units
# panel (b)
Y_AX = -4.75
H_DIR, H_GH = 0.95, 0.32

fig = plt.figure(figsize=(3.45, 3.4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
ax.set_xlim(-1.35, 10.6)
ax.set_ylim(Y_AX - 0.78, 1.0)
ax.set_aspect('auto')
ax.axis('off')

LW_FIB, LW_LANE, LW_GH = 2.6, 1.5, 1.1


def hook(x0, y_from, y_to, x_dir, r=R_HOOK):
    """Quarter-turns from a lane at y_from down to y_to, turning at x0.
    x_dir = -1: the path comes from the left, turns down and leaves to the left
    (a reflection). x_dir = +1: comes from the left, turns down and continues
    to the right (branching)."""
    t = np.linspace(0, np.pi / 2, 12)
    ry = (y_from - y_to) / 2.0
    if x_dir < 0:
        # in from the left at y_from, turn down at x0, out to the left at y_to
        xa = x0 - r + r * np.sin(t)
        ya = y_from - ry + ry * np.cos(t)
        xb = x0 - r + r * np.cos(t)
        yb = y_to + ry - ry * np.sin(t)
        return np.r_[xa, xb], np.r_[ya, yb]
    xa = x0 - r + r * np.sin(t)
    ya = y_from - ry + ry * np.cos(t)
    xb = x0 + r - r * np.cos(t)
    yb = y_to + ry - ry * np.sin(t)
    return np.r_[xa, xb], np.r_[ya, yb]


def arrow(x0, x1, y, color, lw, ms=7):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle='-|>', mutation_scale=ms,
                                 color=color, lw=lw, shrinkA=0, shrinkB=0, zorder=3))


# ---------------- (a) fibre, gratings, launched code ------------------------
ax.plot([X0 + 0.15, XEND], [Y_FIB, Y_FIB], color='0.62', lw=LW_FIB, solid_capstyle='butt', zorder=1)
ax.plot([XEND, XEND + 0.55], [Y_FIB, Y_FIB], color='0.62', lw=LW_FIB, ls=(0, (1.2, 1.4)), zorder=1)
for x, g, col in GR:
    for dx in (-0.075, 0.0, 0.075):
        ax.plot([x + dx, x + dx], [Y_FIB - 0.18, Y_FIB + 0.18], color=col, lw=1.6, solid_capstyle='butt', zorder=2)
    ax.text(x, Y_FIB + 0.3, r'$%s$, $R_%s$' % (g, g), color=col, ha='center', va='bottom', fontsize=7)
ax.text(X0 + 0.15, 0.78, 'upstream', color=GREY, ha='left', va='center', fontsize=6.2)
arrow(X0 + 1.55, X0 + 2.15, 0.78, GREY, 0.7, ms=5)
ax.text(X0 + 2.3, 0.78, 'downstream', color=GREY, ha='left', va='center', fontsize=6.2)
ax.text(XEND + 0.55, 0.78, r'$\tau_k = 2 n_g z_k / c$', color=GREY, ha='right', va='center', fontsize=6.6)

# launched code: runs past every grating
ax.plot([X0 + 0.15, TA], [Y_IN, Y_IN], color='0.5', lw=LW_LANE, zorder=2)
arrow(0.1, 0.55, Y_IN, '0.5', LW_LANE)
arrow(TA - 1.25, TA - 0.75, Y_IN, '0.5', LW_LANE)
ax.text(X0 + 0.05, Y_IN, 'in', color=GREY, ha='right', va='center', fontsize=6.6)

# direct returns: hook at the grating, then straight back to the photodiode
for (x, g, col), y_lane, xa in zip(GR, (Y_RB, Y_RC, Y_RA), (0.35, 0.9, 2.6)):
    hx, hy = hook(x, Y_IN, y_lane, -1)
    ax.plot(hx, hy, color=col, lw=LW_LANE, zorder=3)
    ax.plot([X0 + 0.15, x - R_HOOK], [y_lane, y_lane], color=col, lw=LW_LANE, zorder=2)
    arrow(xa + 0.5, xa, y_lane, col, LW_LANE)
    ax.text(x + 0.1, y_lane - 0.02, r'$\times R_%s$' % g, color=col, ha='left', va='center', fontsize=6.6)

# the ghost (a,b,c): leaves the return of a at b, runs to c, is reflected, returns
gx, gy = hook(TB, Y_RA, Y_G1, +1)
ax.plot(gx, gy, color=PUR, lw=LW_GH, zorder=3)
ax.plot([TB + R_HOOK, TC - R_HOOK], [Y_G1, Y_G1], color=PUR, lw=LW_GH, zorder=2)
arrow(1.5, 1.95, Y_G1, PUR, LW_GH, ms=6)
gx, gy = hook(TC, Y_G1, Y_G2, -1)
ax.plot(gx, gy, color=PUR, lw=LW_GH, zorder=3)
ax.plot([X0 + 0.15, TC - R_HOOK], [Y_G2, Y_G2], color=PUR, lw=LW_GH, zorder=2)
arrow(0.75, 0.3, Y_G2, PUR, LW_GH, ms=6)
ax.text(TB - 0.12, Y_G1 + 0.02, r'$\times R_b$', color=FS.ORANGE, ha='right', va='center', fontsize=6.6)
ax.text(TC + 0.12, Y_G1 + 0.02, r'$\times R_c$', color=FS.GREEN, ha='left', va='center', fontsize=6.6)
ax.text(TC + 0.12, Y_G2 + 0.02, r'path $(a,b,c)$: $\tau_a{-}\tau_b{+}\tau_c$, power $\propto R_aR_bR_c$',
        color=PUR, ha='left', va='center', fontsize=6.4)

# photodiode: everything lands here
ax.plot([X0, X0], [Y_RB + 0.15, Y_G2 - 0.15], color='0.3', lw=1.6, solid_capstyle='butt', zorder=4)
ax.text(X0 - 0.16, (Y_RB + Y_G2) / 2, 'PD', color=GREY, ha='center', va='center', rotation=90, fontsize=6.6)

# ---------------- (b) delay axis, aligned with the fibre --------------------
for x, g, col in GR:
    # leaders: a stub under the grating and a line from below the lanes to the bar,
    # so that they do not run along the vertical part of the hooks
    ax.plot([x, x], [Y_FIB - 0.22, Y_IN + 0.12], color=col, lw=0.7, ls=(0, (1, 1.6)), alpha=0.8, zorder=0)
    ax.plot([x, x], [Y_G2 - 0.3, Y_AX + H_DIR + 0.42], color=col, lw=0.7, ls=(0, (1, 1.6)), alpha=0.8, zorder=0)
arrow(X0 + 0.15, 10.35, Y_AX, '0.3', 0.9, ms=7)
ax.text(10.42, Y_AX, r'$\tau$', color='0.25', ha='left', va='center', fontsize=7)
ax.text(X0 + 0.15, Y_AX + H_DIR + 0.55, 'arrivals at the PD, heights not to scale',
        color=GREY, ha='left', va='center', fontsize=6.2)
for x, g, col in GR:
    ax.plot([x, x], [Y_AX, Y_AX + H_DIR], color=col, lw=3.0, solid_capstyle='butt', zorder=3)
    ax.text(x, Y_AX - 0.14, r'$\tau_%s$' % g, color=col, ha='center', va='top', fontsize=7)
    ax.text(x, Y_AX + H_DIR + 0.05, r'$R_%s$' % g, color=col, ha='center', va='bottom', fontsize=7)
for x, lab_t, lab_r, n in GHOSTS:
    if n == 1:
        ax.plot([x, x], [Y_AX, Y_AX + H_GH], color=PUR, lw=2.2, solid_capstyle='butt', zorder=3)
        ax.text(x, Y_AX - 0.14, lab_t, color=PUR, ha='center', va='top', fontsize=6.4)
        ax.text(x, Y_AX + H_GH + 0.05, lab_r, color=PUR, ha='center', va='bottom', fontsize=6.6)
    else:
        for dx in (-0.09, 0.09):
            ax.plot([x + dx, x + dx], [Y_AX, Y_AX + H_GH], color=PUR, lw=2.2, solid_capstyle='butt', zorder=3)
        ax.text(x, Y_AX - 0.36, lab_t, color=PUR, ha='center', va='top', fontsize=6.4)
        ax.text(x, Y_AX + H_GH + 0.28, lab_r, color=PUR, ha='center', va='bottom', fontsize=6.6)

# panel letters in the house convention
for y, let in ((0.86, 'a'), (Y_AX + H_DIR + 0.55, 'b')):
    ax.text(X0 - 0.55, y, let, fontsize=9, fontweight='bold', ha='left', va='center')

os.makedirs('figs', exist_ok=True)
fig.savefig('figs/fig_s50_ghostpath.pdf')
fig.savefig('figs/fig_s50_ghostpath.png', dpi=220)
print('saved figs/fig_s50_ghostpath.pdf')
