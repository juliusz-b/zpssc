"""s53_ruler.py - gratings on a Golomb ruler, in the style of the ghost figure.

Replaces the TikZ drawing fig_ruler.tex.

  (a) The ruler {0, 1, 4, 6}: four gratings on its marks, and the six
      pairwise differences drawn as brackets, each value occurring once.
  (b) The marks as delay bins of a periodic code with N = 13. Tall bars are
      the gratings, short bars every third-order ghost at its delay
      tau_i - tau_j + tau_l, in the greyed colour of its first grating,
      side by side where several share a bin. None lands on a grating.
  (c) The same marks with N = 11: the delays 11 and 12 fold back to bins
      0 and 1, and two ghosts land on gratings.
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

MARKS = [0, 1, 4, 6]
COL = [FS.ORANGE, FS.GREEN, FS.VERM, FS.BLUE]      # gratings 1..4 by position
GREY = '0.45'


def tint(col, w):
    c = np.array(to_rgb(col)); g = np.array([0.78, 0.78, 0.78])
    return tuple(w * c + (1 - w) * g)


GH = [tint(c, 0.45) for c in COL]


def ghosts():
    """[(delay, index of the first grating), ...] for every third-order path."""
    out = []
    K = len(MARKS)
    for i in range(K):
        for j in range(K):
            for l in range(K):
                if MARKS[j] < MARKS[i] and MARKS[j] < MARKS[l]:
                    out.append((MARKS[i] - MARKS[j] + MARKS[l], i))
    return out


fig = plt.figure(figsize=(3.45, 2.8))
gs = fig.add_gridspec(3, 1, height_ratios=[1.25, 1.0, 1.0], hspace=0.45,
                      left=0.03, right=0.9, bottom=0.09, top=0.93)
ax, bx, cx = (fig.add_subplot(gs[i, 0]) for i in range(3))

# ---------------- (a) the ruler and its differences -------------------------
ax.plot([0, 6], [0, 0], color='0.3', lw=0.9, zorder=1)
for m in range(7):
    ax.plot([m, m], [0, -0.08], color='0.3', lw=0.6)
    ax.text(m, -0.14, str(m), ha='center', va='top', fontsize=6.2, color='0.3')
for g, m in enumerate(MARKS):
    ax.plot([m, m], [0, 0.3], color=COL[g], lw=2.6, solid_capstyle='butt', zorder=3)
pairs = sorted(((MARKS[q] - MARKS[p], p, q) for p in range(4) for q in range(p + 1, 4)))
for n, (d, p, q) in enumerate(pairs):
    h = 0.45 + 0.19 * n
    x0, x1 = MARKS[p], MARKS[q]
    ax.plot([x0, x0, x1, x1], [0.34, h, h, 0.34], color='0.55', lw=0.6, zorder=2)
    ax.text(0.5 * (x0 + x1), h, str(d), ha='center', va='center', fontsize=6.2, color='0.25',
            bbox=dict(facecolor='white', edgecolor='none', pad=0.8), zorder=4)
ax.text(6.35, 0.15, 'delay bin', color=GREY, ha='left', va='center', fontsize=6.0)
ax.set_xlim(-0.4, 7.6)
ax.set_ylim(-0.35, 1.55)
ax.axis('off')
ax.text(0.035, 1.0, r'ruler $\{0,1,4,6\}$: every difference occurs once', transform=ax.transAxes,
        ha='left', va='bottom', fontsize=6.4, color=GREY)


# ---------------- (b), (c) delay bins of a periodic code --------------------
def bins_panel(a_, N, note):
    a_.plot([0, N - 1], [0, 0], color='0.3', lw=0.9, zorder=1)
    a_.plot([N - 1, N - 0.3], [0, 0], color='0.3', lw=0.9, ls=(0, (1.2, 1.4)), zorder=1)
    a_.add_patch(FancyArrowPatch((N - 0.55, 0.0), (N - 0.55, -0.55), connectionstyle='arc3,rad=-0.6',
                                 arrowstyle='-|>', mutation_scale=6, color='0.4', lw=0.7, zorder=1))
    a_.text(N - 0.35, -0.5, 'to 0', color=GREY, ha='left', va='center', fontsize=5.8)
    for m in range(N):
        a_.plot([m, m], [0, -0.08], color='0.3', lw=0.6)
        a_.text(m, -0.14, str(m), ha='center', va='top', fontsize=6.0, color='0.3')
    for g, m in enumerate(MARKS):
        a_.plot([m, m], [0, 1.0], color=COL[g], lw=2.6, solid_capstyle='butt', zorder=3)
    per_bin = {}
    for d, g in ghosts():
        per_bin.setdefault(d % N, []).append((d, g))
    for b_, lst in per_bin.items():
        n = len(lst)
        for k, (d, g) in enumerate(lst):
            off = (k - 0.5 * (n - 1)) * 0.16
            if d >= N:
                off += 0.22                             # beside the grating bar it lands on
            a_.plot([b_ + off, b_ + off], [0, 0.38], color=GH[g], lw=2.0, solid_capstyle='butt', zorder=2)
            if d >= N:                                   # folded onto a grating
                a_.plot(b_ + off, 0.52, 'x', color=FS.VERM, ms=4.5, mew=1.1, zorder=5)
    a_.set_xlim(-0.4, N + 1.4)
    a_.set_ylim(-0.6, 1.25)
    a_.axis('off')
    a_.text(0.035, 1.0, note, transform=a_.transAxes, ha='left', va='bottom', fontsize=6.4, color=GREY)


bins_panel(bx, 13, r'delay bins of a code with $N=13$: no ghost on a grating')
bins_panel(cx, 11, r'$N=11$: two ghosts fold back onto gratings')
cx.text(1.9, 0.75, r'$11\to0$, $12\to1$', color=FS.VERM, ha='left', va='center', fontsize=6.0)

for a_, let in ((ax, 'a'), (bx, 'b'), (cx, 'c')):
    a_.text(-0.01, 1.02, let, transform=a_.transAxes, fontsize=9, fontweight='bold', ha='right', va='bottom')

os.makedirs('figs', exist_ok=True)
fig.savefig('figs/fig_s53_ruler.pdf')
fig.savefig('figs/fig_s53_ruler.png', dpi=220)
print('saved figs/fig_s53_ruler.pdf')
print('ghost delays:', sorted(d for d, _ in ghosts()))
