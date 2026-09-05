"""s53_ruler.py - why gratings on a Golomb ruler never receive a ghost.

Replaces the TikZ drawing fig_ruler.tex. Blue squares mark gratings, grey
dots mark third-order paths at delay tau_i - tau_j + tau_l, and orange
crosses mark paths arriving at occupied bins. One marker represents one
path, with coincident paths stacked vertically at the exact delay bin.
The delay axis is in chips of the code, one bin per chip.

  (a) Four gratings at uniform spacing {0, 1, 2, 3}: every ghost delay is
      again a grid point, and from the third grating on the ghosts land on
      gratings (crosses).
  (b) The same four gratings on the Golomb ruler {0, 1, 4, 6}. The six
      pairwise differences, drawn as brackets, each occur once, so no ghost
      lands on a grating.
  (c) The ruler with a code of N = 7 chips: every delay beyond 6 folds back
      and six ghosts land on the gratings at 0, 1 and 4. Both N are real
      m-sequence lengths.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import to_rgb
import figstyle as FS

FS.apply()

UNIFORM = [0, 1, 2, 3]
RULER = [0, 1, 4, 6]
COL = [FS.ORANGE, FS.GREEN, FS.VERM, FS.BLUE]      # gratings 1..4 by position
GREY = '0.45'


def tint(col, w):
    c = np.array(to_rgb(col)); g = np.array([0.78, 0.78, 0.78])
    return tuple(w * c + (1 - w) * g)


GH = [tint(c, 0.45) for c in COL]


def ghosts(marks):
    """[(delay, index of the first grating), ...] for every third-order path."""
    out = []
    K = len(marks)
    for i in range(K):
        for j in range(K):
            for l in range(K):
                if marks[j] < marks[i] and marks[j] < marks[l]:
                    out.append((marks[i] - marks[j] + marks[l], i))
    return out


# Separate physical grating positions from the count of paths at each delay.
from matplotlib.lines import Line2D

fig = plt.figure(figsize=(3.45, 3.45))
gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.75, 1.0], hspace=0.54,
                      left=0.06, right=0.98, bottom=0.09, top=0.84)
ax, bx, cx = (fig.add_subplot(gs[i, 0]) for i in range(3))
legend = [Line2D([], [], marker='s', ls='none', color=FS.BLUE, ms=4,
                 label='Grating'),
          Line2D([], [], marker='o', ls='none', color='0.50', ms=3.2,
                 label='Ghost path'),
          Line2D([], [], marker='x', ls='none', color=FS.VERM, ms=4,
                 label='Collision')]
fig.legend(handles=legend, loc='upper center', bbox_to_anchor=(0.52, 1.005),
           ncol=3, frameon=False, fontsize=6.2, handlelength=0.8,
           handletextpad=0.4, columnspacing=0.9)


def bins_panel(axis, marks, N, title, letter, brackets=False):
    axis.set_xlim(-0.7, 14.7)
    axis.set_ylim(-1.42, 1.77 if brackets else 0.28)
    axis.axis('off')
    axis.text(-0.035, 1.10, letter, transform=axis.transAxes,
              fontweight='bold', fontsize=9, va='bottom')
    axis.text(0.04, 1.10, title, transform=axis.transAxes,
              fontsize=6.9, fontweight='bold', va='bottom', color='0.22')
    axis.plot([-0.35, N - 0.65], [0, 0], color='0.35', lw=0.85)
    for m in range(N):
        axis.plot([m, m], [0, -1.12], color='0.90', lw=0.45, zorder=0)
        axis.text(m, -1.22, str(m), ha='center', va='top', fontsize=6.0,
                  color='0.30')
    for m in marks:
        axis.plot(m, 0, 's', color=FS.BLUE, ms=4.3, zorder=4)
    per_bin = {}
    for delay, first in ghosts(marks):
        per_bin.setdefault(delay % N, []).append((delay, first))
    hits = 0
    for delay, paths_here in sorted(per_bin.items()):
        collision = delay in marks
        if collision:
            hits += len(paths_here)
        for j, _ in enumerate(paths_here):
            axis.plot(delay, -0.25 - 0.17 * j,
                      marker='x' if collision else 'o', ls='none',
                      ms=3.3 if collision else 2.8, mew=0.95,
                      color=FS.VERM if collision else '0.50', zorder=3)
    if brackets:
        pairs = sorted((marks[j] - marks[i], marks[i], marks[j])
                       for i in range(len(marks)) for j in range(i + 1, len(marks)))
        for k, (distance, x0, x1) in enumerate(pairs):
            y = 0.28 + 0.24 * k
            axis.plot([x0, x0, x1, x1], [0.10, y, y, 0.10],
                      color='0.60', lw=0.6, zorder=1)
            axis.text((x0 + x1) / 2, y, str(distance), fontsize=6.0,
                      ha='center', va='center', color='0.25',
                      bbox=dict(facecolor='white', edgecolor='none', pad=0.4))
    score_y = 1.03 if brackets else -0.20
    axis.text(10.2, score_y, '%d of 14 paths collide' % hits,
              fontsize=6.2, color=FS.VERM if hits else FS.BLUE,
              ha='center', va='center')
    if N == 7:
        axis.axvline(6.5, ymin=0.08, ymax=0.88, lw=0.7,
                     color='0.65', ls=(0, (2, 2)))
        axis.text(10.2, -0.59, 'Wrapped delays', ha='center', fontsize=5.8,
                  color='0.40')
        axis.text(10.2, -0.94, r'$7\to0,\quad8\to1,\quad11\to4$',
                  ha='center', fontsize=6.0, color=FS.VERM)
    return hits


h_a = bins_panel(ax, UNIFORM, 15, 'Uniform spacing, $N=15$', 'a')
h_b = bins_panel(bx, RULER, 15, 'Golomb ruler, $N=15$', 'b', brackets=True)
h_c = bins_panel(cx, RULER, 7, 'Same ruler, $N=7$', 'c')
fig.text(0.52, 0.022, 'Delay bin (chips)', ha='center', fontsize=7.0)
print('ghosts on gratings: uniform %d, ruler %d, ruler with N=7 %d' % (h_a, h_b, h_c))
os.makedirs('figs', exist_ok=True)
fig.savefig('figs/fig_s53_ruler.pdf')
fig.savefig('figs/fig_s53_ruler.png', dpi=300)
print('saved figs/fig_s53_ruler.pdf')
