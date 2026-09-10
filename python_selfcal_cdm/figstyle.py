"""figstyle.py - one visual language for every figure in the paper.

Figures are drawn at their final printed width (7.1 in for a full-width IEEE
figure), so a 7 pt label on screen is a 7 pt label on paper. That holds only
while the PDF width matches the width main.tex includes it at, which is worth
checking after any layout change: pdfinfo gives the file width, and the ratio
against the includegraphics width multiplies every font size on the figure. Palette is
Paul Tol's vibrant set, colour-blind safe, with fixed roles (C_*) shared by all figures. Since 31.08.2026 the look follows Origin: boxed
axes, ticks inside, sans-serif labels, framed legend. Panels lettered in bold
lowercase.
"""
import matplotlib.pyplot as plt

# Paul Tol "vibrant" palette (colour-blind safe, high contrast). The old
# Okabe-Ito names are kept so that every script switches at once.
BLUE = '#0077BB'
SKY = '#33BBEE'
ORANGE = '#EE7733'
VERM = '#CC3311'
GREEN = '#009988'
PURPLE = '#EE3377'
YELLOW = '#EE7733'
GREY = '#4D4D4D'
LGREY = '#BBBBBB'

# Roles: the same quantity gets the same colour in every figure of the paper.
C_TRUE = '#9A9A9A'      # isolated (true) spectrum R_k, drawn thick
C_THEORY = '#222222'    # analytical rule or bound
C_DIRECT = BLUE         # direct return A_k, direct paths, gratings on the delay axis
C_MEAS = VERM           # measured or uncorrected, initial array, narrow first, uniform spacing, 4 m, R = 10 %
C_CORR = GREEN          # corrected, deshadowed
C_GOOD = BLUE           # randomized spacing, wide first, 40 m, R = 1 %, designed array
C_GHOST = ORANGE        # ghosts and their paths
C_LEAK = '#777777'      # code leakage
C_SHADOW = PURPLE       # spectral shadowing
C_STEP = (ORANGE, SKY)  # intermediate correction steps
REFS = {0: VERM, 1: ORANGE, 2: BLUE, 3: GREEN}   # number of reference gratings


def apply(base=7.0):
    """Styl jak z Origina (decyzja kierownika, 31.08.2026): ramka osi z
    czterech stron, znaczniki glowne i pomocnicze do wewnatrz, krój
    bezszeryfowy, pogrubione opisy osi, grubsze linie, legenda w ramce.
    Poprzedni styl (Times, otwarte osie) zostal w apply_open()."""
    plt.rcParams.update({
        'font.size': base,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
        'mathtext.fontset': 'dejavusans',
        'axes.titlesize': base + 1,
        'axes.labelsize': base + 1,
        'axes.labelweight': 'bold',
        'xtick.labelsize': base,
        'ytick.labelsize': base,
        'legend.fontsize': base - 1.0,
        'legend.handlelength': 1.8,
        'axes.grid': False,
        'axes.linewidth': 1.0,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.top': True,
        'ytick.right': True,
        'xtick.minor.visible': True,
        'ytick.minor.visible': True,
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'xtick.minor.size': 2,
        'ytick.minor.size': 2,
        'xtick.major.width': 1.0,
        'ytick.major.width': 1.0,
        'xtick.minor.width': 0.7,
        'ytick.minor.width': 0.7,
        'lines.linewidth': 1.4,
        'lines.markersize': 4.0,
        'axes.spines.top': True,
        'axes.spines.right': True,
        'legend.frameon': True,
        'legend.fancybox': False,
        'legend.edgecolor': 'black',
        'legend.framealpha': 1.0,
        'figure.dpi': 200,
        'savefig.dpi': 300,
    })


def apply_open(base=7.0):
    """Poprzedni styl: Times jak tekst IEEEtran, cienkie otwarte osie,
    legenda bez ramki. Dla figur, ktorych jeszcze nie przerobiono."""
    plt.rcParams.update({
        'font.size': base,
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


def letter(ax, s):
    """Panel letter in the convention of Fig. 4: bold lowercase, above the top-left corner."""
    ax.text(0.02, 1.06, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='bottom')


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


# ---------------------------------------------------------------------------
# Axis labels start with a capital letter (decision of 6.09.2026). Applied at
# import, so every script that uses this module gets it without edits. Labels
# that start with math ($...$) or a symbol are left alone.
# ---------------------------------------------------------------------------
import matplotlib.axes as _mpl_axes


def _capitalize_label(s):
    if isinstance(s, str) and s and s[0].islower():
        return s[0].upper() + s[1:]
    return s


_orig_set_xlabel = _mpl_axes.Axes.set_xlabel
_orig_set_ylabel = _mpl_axes.Axes.set_ylabel


def _set_xlabel(self, xlabel, *a, **k):
    return _orig_set_xlabel(self, _capitalize_label(xlabel), *a, **k)


def _set_ylabel(self, ylabel, *a, **k):
    return _orig_set_ylabel(self, _capitalize_label(ylabel), *a, **k)


_mpl_axes.Axes.set_xlabel = _set_xlabel
_mpl_axes.Axes.set_ylabel = _set_ylabel
