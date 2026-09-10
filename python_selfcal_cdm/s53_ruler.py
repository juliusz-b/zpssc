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

FS.apply(base=6.5)
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


fig = plt.figure(figsize=(3.45, 2.25))
axes = [fig.add_axes([0.075, bottom, 0.905, 0.185])
        for bottom in (0.66, 0.385, 0.11)]
legend = [Line2D([], [], marker='s', ls='none', color=FS.BLUE, ms=3.8,
                 label='Grating'),
          Line2D([], [], marker='o', ls='none', color='0.50', ms=2.5,
                 label='Ghost path'),
          Line2D([], [], marker='x', ls='none', color=FS.VERM, ms=3.4,
                 label='Collision')]
fig.legend(handles=legend, loc='upper center', bbox_to_anchor=(0.53, 1.005),
           ncol=3, frameon=False, fontsize=6.1, handlelength=0.8,
           handletextpad=0.4, columnspacing=0.9)


def bins_panel(axis, marks, period, title, letter):
    axis.set_xlim(-0.35, 14.35)
    axis.set_ylim(-1.40, 0.20)
    axis.axis('off')
    paths = ghosts(marks)
    counts = Counter(delay % period for delay, _ in paths)
    hits = sum(count for delay, count in counts.items() if delay in marks)
    axis.text(-0.067, 1.07, letter, transform=axis.transAxes,
              fontweight='bold', fontsize=8, va='bottom')
    axis.text(0, 1.07, title, transform=axis.transAxes,
              fontsize=6.3, va='bottom', color='0.20')
    axis.text(1, 1.07, f'{hits}/{len(paths)} collide', transform=axis.transAxes,
              fontsize=6.0, va='bottom', ha='right',
              color=FS.VERM if hits else FS.BLUE)
    axis.plot([-0.2, period - 0.7], [0, 0], color='0.35', lw=0.75)
    for delay in range(period):
        axis.plot([delay, delay], [0, -1.03], color='0.9', lw=0.35, zorder=0)
        axis.text(delay, -1.12, str(delay), ha='center', va='top',
                  fontsize=5.5, color='0.3')
    axis.plot(marks, [0]*len(marks), 's', color=FS.BLUE, ms=3.7, zorder=3)
    for delay, count in sorted(counts.items()):
        collision = delay in marks
        for j in range(count):
            axis.plot(delay, -0.23-0.15*j, marker='x' if collision else 'o',
                      ls='none', ms=2.8 if collision else 2.3, mew=0.8,
                      color=FS.VERM if collision else '0.5', zorder=3)
    return dict(marks=marks, period=period, paths=paths,
                counts=dict(sorted(counts.items())), collisions=hits)


report = [bins_panel(axes[0], UNIFORM, 15, 'Uniform spacing, $N=15$', 'a'),
          bins_panel(axes[1], RULER, 15, 'Golomb ruler, $N=15$', 'b'),
          bins_panel(axes[2], RULER, 7, 'Same ruler, $N=7$', 'c')]

c = axes[2]
c.plot([6.5, 6.5], [0.12, -1.03], color='0.60', lw=0.65, ls=(0,(2,2)))
c.text(10.6, -0.15, 'Wrapped delays', fontsize=5.6, color='0.30',
       ha='center', va='top')
c.text(10.6, -0.60, r'$7\to0,\quad8\to1,\quad11\to4$',
       fontsize=6.0, color=FS.VERM, ha='center', va='center')
fig.text(0.53, 0.023, 'Delay bin [chips]', ha='center', fontsize=6.5)

assert [r['collisions'] for r in report] == [4, 0, 6]
assert all(len(r['paths']) == 14 for r in report)
OUT.mkdir(exist_ok=True)
fig.savefig(OUT/'fig_s53_ruler.pdf')
fig.savefig(OUT/'fig_s53_ruler.png', dpi=300)
(OUT/'s53_ruler_results.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print('Collisions: 4, 0, 6 of 14 paths.')
