"""figstyle.py - one visual language for every figure in the paper.

Figures are drawn at their final printed width (7.1 in for a full-width IEEE
figure), so a 7 pt label on screen is a 7 pt label on paper. Palette is
Okabe-Ito, colour-blind safe. Spines are thin and open (top and right off),
legends frameless, panels lettered in bold lowercase.
"""
import matplotlib.pyplot as plt

# Okabe-Ito
BLUE = '#0072B2'
SKY = '#56B4E9'
ORANGE = '#E69F00'
VERM = '#D55E00'
GREEN = '#009E73'
PURPLE = '#CC79A7'
YELLOW = '#F0E442'
GREY = '#4D4D4D'
LGREY = '#AAAAAA'


def apply(base=7.0):
    plt.rcParams.update({
        'font.size': base,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
        'mathtext.fontset': 'dejavusans',
        'axes.titlesize': base,
        'axes.labelsize': base,
        'xtick.labelsize': base - 0.5,
        'ytick.labelsize': base - 0.5,
        'legend.fontsize': base - 0.5,
        'axes.linewidth': 0.6,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 2.2,
        'ytick.major.size': 2.2,
        'lines.linewidth': 1.0,
        'lines.markersize': 3.2,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'legend.frameon': False,
        'figure.dpi': 200,
        'savefig.dpi': 300,
    })


def panel(fig, x, y, letter, size=8.5):
    """Bold lowercase panel letter at figure coordinates."""
    fig.text(x, y, letter, fontsize=size, fontweight='bold', va='top')


def despine_all(ax):
    for s in ('top', 'right', 'left', 'bottom'):
        ax.spines[s].set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
