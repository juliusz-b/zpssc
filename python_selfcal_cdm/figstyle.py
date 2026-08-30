"""figstyle.py - one visual language for every figure in the paper.

Figures are drawn at their final printed width (7.1 in for a full-width IEEE
figure), so a 7 pt label on screen is a 7 pt label on paper. That holds only
while the PDF width matches the width main.tex includes it at, which is worth
checking after any layout change: pdfinfo gives the file width, and the ratio
against the includegraphics width multiplies every font size on the figure. Palette is
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
        # Times, jak tekst artykulu skladany IEEEtranem. Figury przestaja
        # wygladac jak wklejone z innego dokumentu, a wzory w podpisach osi
        # maja ten sam krój co wzory w tekscie. STIX to wolny odpowiednik
        # Timesa z pelnym zestawem matematycznym.
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'STIXGeneral', 'DejaVu Serif'],
        'mathtext.fontset': 'stix',
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


def dim_gap(ax, x1, x2, y, label, color='k', tail=None, side='right',
            fontsize=6.2, lw=0.8, head=5.5, pad=None):
    """Wymiaruj odstep za maly, zeby zmiescic w nim groty strzalek.

    Tak sie to robi na rysunku technicznym: skoro strzalki nie mieszcza sie
    miedzy liniami wymiarowymi, wychodza na zewnatrz i pokazuja do srodka, a
    liczba przenosi sie obok. Same granice odstepu musza byc juz zaznaczone,
    zwykle pionowymi liniami krzywych, ktore sie porownuje.

    tail to dlugosc strzalki na zewnatrz, domyslnie dziesiata czesc osi.
    side mowi, po ktorej stronie zostalo miejsce na liczbe.
    """
    lo, hi = (x1, x2) if x1 <= x2 else (x2, x1)
    span = ax.get_xlim()[1] - ax.get_xlim()[0]
    if tail is None:
        tail = 0.10 * span
    if pad is None:
        pad = 0.35 * tail
    ax.plot([lo, hi], [y, y], color=color, lw=lw, solid_capstyle='butt',
            zorder=6)
    for x, sgn in ((lo, -1.0), (hi, 1.0)):
        ax.annotate('', xy=(x, y), xytext=(x + sgn * tail, y),
                    arrowprops=dict(arrowstyle='-|>', color=color, lw=lw,
                                    mutation_scale=head, shrinkA=0,
                                    shrinkB=0), zorder=6)
    if side == 'right':
        ax.text(hi + tail + pad, y, label, ha='left', va='center',
                fontsize=fontsize, color=color, zorder=6)
    else:
        ax.text(lo - tail - pad, y, label, ha='right', va='center',
                fontsize=fontsize, color=color, zorder=6)
