"""s53_ruler.py - why gratings on a Golomb ruler never receive a ghost.

Replaces the TikZ drawing fig_ruler.tex. Same vocabulary as the ghost
figure: tall bars are gratings, short bars are third-order ghosts at their
delay tau_i - tau_j + tau_l, in the greyed colour of the grating of their
first reflection. The delay axis is in chips of the code, one bin per chip.

  (a) Four gratings at uniform spacing {0, 1, 2, 3}: every ghost delay is
      again a grid point, and from the third grating on the ghosts land on
      gratings (crosses).
  (b) The same four gratings on the Golomb ruler {0, 1, 4, 6}. The six
      pairwise differences, drawn as brackets, each occur once, so no ghost
      lands on a grating.
  (c) The ruler with a code of N = 11 chips: the delays 11 and 12 fold back
      to bins 0 and 1, and two ghosts land on gratings.
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


fig = plt.figure(figsize=(3.45, 3.1))
gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.55, 1.0], hspace=0.5,
                      left=0.03, right=0.9, bottom=0.08, top=0.94)
ax, bx, cx = (fig.add_subplot(gs[i, 0]) for i in range(3))
XMAX = 13


def bins_panel(a_, marks, N, note, brackets=False):
    """Delay axis of N bins with the gratings and every ghost."""
    a_.plot([0, N - 1], [0, 0], color='0.3', lw=0.9, zorder=1)
    if N <= XMAX:
        a_.plot([N - 1, N - 0.3], [0, 0], color='0.3', lw=0.9, ls=(0, (1.2, 1.4)), zorder=1)
        a_.add_patch(FancyArrowPatch((N - 0.55, 0.0), (N - 0.55, -0.55), connectionstyle='arc3,rad=-0.6',
                                     arrowstyle='-|>', mutation_scale=6, color='0.4', lw=0.7, zorder=1))
        a_.text(N - 0.35, -0.5, 'to 0', color=GREY, ha='left', va='center', fontsize=5.8)
    for m in range(N):
        a_.plot([m, m], [0, -0.08], color='0.3', lw=0.6)
        a_.text(m, -0.14, str(m), ha='center', va='top', fontsize=6.0, color='0.3')
    for g, m in enumerate(marks):
        a_.plot([m, m], [0, 1.0], color=COL[g], lw=2.6, solid_capstyle='butt', zorder=3)
    per_bin = {}
    for d, g in ghosts(marks):
        per_bin.setdefault(d % N, []).append((d, g))
    hits = 0
    for b_, lst in per_bin.items():
        n = len(lst)
        for k, (d, g) in enumerate(lst):
            off = (k - 0.5 * (n - 1)) * 0.16
            on_grating = b_ in marks
            if on_grating:
                off += 0.24                              # beside the grating bar it lands on
            a_.plot([b_ + off, b_ + off], [0, 0.38], color=GH[g], lw=2.0, solid_capstyle='butt', zorder=2)
            if on_grating:
                hits += 1
        if b_ in marks:
            a_.plot(b_, 1.14, 'x', color=FS.VERM, ms=5, mew=1.2, zorder=5, clip_on=False)
    if brackets:
        pairs = sorted(((marks[q] - marks[p], p, q) for p in range(4) for q in range(p + 1, 4)))
        for n, (d, p, q) in enumerate(pairs):
            h = 1.25 + 0.24 * n
            x0, x1 = marks[p], marks[q]
            a_.plot([x0, x0, x1, x1], [1.06, h, h, 1.06], color='0.55', lw=0.6, zorder=2)
            a_.text(0.5 * (x0 + x1), h, str(d), ha='center', va='center', fontsize=6.0, color='0.25',
                    bbox=dict(facecolor='white', edgecolor='none', pad=0.6), zorder=4)
    a_.set_xlim(-0.4, XMAX + 1.4)
    a_.set_ylim(-0.6, 2.75 if brackets else 1.35)
    a_.axis('off')
    a_.text(0.035, 1.0, note, transform=a_.transAxes, ha='left', va='bottom', fontsize=6.4, color=GREY)
    return hits


h_a = bins_panel(ax, UNIFORM, 13, r'uniform spacing $\{0,1,2,3\}$: ghosts land on gratings')
h_b = bins_panel(bx, RULER, 13, r'Golomb ruler $\{0,1,4,6\}$: every difference once, no ghost on a grating', brackets=True)
h_c = bins_panel(cx, RULER, 11, r'the ruler with $N=11$: two ghosts fold back onto gratings')
cx.text(1.9, 0.75, r'$11\to0$, $12\to1$', color=FS.VERM, ha='left', va='center', fontsize=6.0)
cx.text(XMAX + 1.4, -0.14, 'delay (chips)', color=GREY, ha='right', va='top', fontsize=6.0)
print('ghosts on gratings: uniform %d, ruler %d, ruler with N=11 %d' % (h_a, h_b, h_c))

for a_, let in ((ax, 'a'), (bx, 'b'), (cx, 'c')):
    a_.text(-0.01, 1.02, let, transform=a_.transAxes, fontsize=9, fontweight='bold', ha='right', va='bottom')

os.makedirs('figs', exist_ok=True)
fig.savefig('figs/fig_s53_ruler.pdf')
fig.savefig('figs/fig_s53_ruler.png', dpi=220)
print('saved figs/fig_s53_ruler.pdf')
