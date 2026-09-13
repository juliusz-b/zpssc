"""Exact third-order ghost collisions for periodic codes.

One marker per directed path. Coincident paths are stacked at the same bin.
The three compact rows show direct returns, ghost multiplicity and wrapping.
"""
from pathlib import Path
from collections import Counter
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import figstyle as FS

FS.apply()
plt.rcParams.update({'pdf.fonttype': 42})
UNIFORM = [0, 1, 2, 3]
RULER = [0, 1, 4, 6]
OUT = Path(__file__).resolve().parent / 'figs'

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


from matplotlib.patches import Rectangle

CASES = [(UNIFORM, 15, 'Uniform spacing, $N=15$', 'a'),
         (RULER, 15, 'Golomb ruler, $N=15$', 'b'),
         (RULER, 7, 'Same ruler, $N=7$', 'c')]
ROWS = [1.0, 0.0, -1.0]


def bins_row(axis, y0, marks, period, title, letter):
    """One row of the figure: baseline with the gratings, one marker per ghost path stacked below its bin,
    grating bins shaded so that collisions (crosses inside the shading) stand out."""
    paths = ghosts(marks)
    counts = Counter(delay % period for delay, _ in paths)
    hits = sum(count for delay, count in counts.items() if delay in marks)
    for m in marks:
        axis.add_patch(Rectangle((m - 0.5, y0 - 0.42), 1.0, 0.62, facecolor=FS.BLUE, alpha=0.12,
                                 edgecolor='none', zorder=0))
    axis.plot([-0.3, period - 0.7], [y0, y0], color='0.25', lw=1.0)
    axis.plot(marks, [y0] * len(marks), 's', color=FS.BLUE, ms=5.0, zorder=3)
    for delay, count in sorted(counts.items()):
        collision = delay in marks
        for j in range(count):
            axis.plot(delay, y0 - 0.13 - 0.085 * j, marker='x' if collision else 'o', ls='none',
                      ms=3.6 if collision else 2.8, mew=1.0, mfc=(FS.VERM if collision else 'white'),
                      color=FS.VERM if collision else FS.C_GHOST, zorder=3)
    axis.text(-0.4, y0 + 0.3, title, fontsize=6.5, va='center', ha='left')
    axis.text(14.5, y0 + 0.3, f'{hits}/{len(paths)} collide', fontsize=6.5, va='center', ha='right',
              color=FS.VERM if hits else FS.BLUE)
    axis.text(-1.35, y0 + 0.3, letter, fontsize=8, fontweight='bold', va='center', ha='center', clip_on=False)
    if period < 15:
        axis.plot([period - 0.5, period - 0.5], [y0 - 0.45, y0 + 0.15], color='0.5', lw=0.8, ls=(0, (2, 2)))
        axis.text(10.5, y0 - 0.2, 'wrapped: $7{\\to}0$, $8{\\to}1$, $11{\\to}4$', fontsize=5.8,
                  color=FS.VERM, ha='center', va='center')
    return dict(marks=marks, period=period, paths=paths,
                counts=dict(sorted(counts.items())), collisions=hits)


fig, ax = plt.subplots(figsize=(3.5, 2.3), layout='constrained')
report = [bins_row(ax, y0, marks, period, title, letter) for y0, (marks, period, title, letter) in zip(ROWS, CASES)]
ax.set_xlim(-0.6, 14.6)
ax.set_ylim(-1.55, 1.45)
ax.set_xticks(range(15))
ax.set_yticks([])
ax.tick_params(axis='y', length=0)
ax.set_xlabel('Delay bin [chips]')
legend = [Line2D([], [], marker='s', ls='none', color=FS.BLUE, ms=5, label='grating'),
          Line2D([], [], marker='o', ls='none', color=FS.C_GHOST, ms=2.8, mfc='white', mew=1.0, label='ghost path'),
          Line2D([], [], marker='x', ls='none', color=FS.VERM, ms=3.6, mew=1.0, label='collision')]
ax.legend(handles=legend, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=6,
          handlelength=1.0, columnspacing=0.9, handletextpad=0.4, frameon=False)

assert [r['collisions'] for r in report] == [4, 0, 6]
assert all(len(r['paths']) == 14 for r in report)
OUT.mkdir(exist_ok=True)
fig.savefig(OUT/'fig_s53_ruler.pdf')
fig.savefig(OUT/'fig_s53_ruler.png', dpi=300)
(OUT/'s53_ruler_results.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print('Collisions: 4, 0, 6 of 14 paths.')
