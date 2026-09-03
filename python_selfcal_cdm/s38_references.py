"""s38_references.py - where the reference gratings go, and why it matters.

The paper asks for stabilized reference gratings and says they belong nearest
the interrogator so that nothing shadows them. That is true of the first one
only. The second looks through the first, the third through two, and every
sensor through all three. The axis correction is therefore built on points
that crossed a different number of gratings than the points it is applied to,
and it does not close over the difference.

Worse, the references have to be spread across the band, because the error
being fitted is a function of band position. Spreading them is exactly what
puts them at the detunings where Law A bites hardest. At the procured 10
percent reflectivity the three references disagree with each other by more
than the error they are there to remove, and the correction ends up leaving
more error than doing nothing at all.

A coupler fixes it. Put the references on a short stub of their own and they
shadow nothing and nothing shadows them, so the fit measures the source error
and only the source error. The cost is 6 dB on the round trip, which a link
budget with 99 dB of margin does not notice.

This is a design rule the paper was missing, and it follows from its own
Law A without any new physics.

Output: figs/fig_s38_references.pdf, full text width.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
import figstyle as FS

FS.apply(7.0)

VERM, ORAN, GREE, PURP = '#D55E00', '#E69F00', '#009E73', '#CC79A7'
BLUE, GREY = '#0072B2', '0.45'

FWHM = 250.0
SIG = FWHM / 2.35482
CA = 4.0 / 3.0 * np.sqrt(2.0 / 3.0)

REFS = np.array([-180.0, 0.0, 180.0])          # band positions, pm
SENS = np.linspace(-200.0, 200.0, 9)
AXIS_PEAK = 6.0                                 # assumed smooth axis error, pm


def law_a(dl, r):
    """Law A: shift of the fitted center caused by one upstream neighbor."""
    return -CA * r * dl * np.exp(-dl ** 2 / (3.0 * SIG ** 2))


def axis_error(nu):
    """The smooth wavelength-axis error the references exist to remove."""
    return AXIS_PEAK * nu / 200.0


def residual(r_ref, inline=True):
    """Sensor error left after the reference fit, in pm."""
    if not inline:
        # references on their own stub: nothing shadows them, they shadow
        # nothing, so the fit sees the axis error alone and cancels it
        return np.zeros_like(SENS)
    read = np.array([axis_error(r) + sum(law_a(REFS[j] - r, r_ref)
                                         for j in range(i))
                     for i, r in enumerate(REFS)])
    coef = np.polyfit(REFS, read, len(REFS) - 1)
    seen = np.array([axis_error(s) + sum(law_a(r - s, r_ref) for r in REFS)
                     for s in SENS])
    return seen - np.polyval(coef, SENS)


fig = plt.figure(figsize=(7.16, 2.30))
axa = fig.add_axes([0.010, 0.06, 0.400, 0.86])
axb = fig.add_axes([0.492, 0.175, 0.222, 0.745])
axc = fig.add_axes([0.792, 0.175, 0.200, 0.745])

# ==========================================================================
# (a) the two places a reference can sit
# ==========================================================================
axa.set_xlim(0, 10.6)
axa.set_ylim(0.30, 5.35)
axa.axis('off')


def source(y):
    axa.add_patch(Rectangle((0.05, y - 0.26), 1.15, 0.52, fc='white',
                            ec='0.35', lw=0.8))
    axa.text(0.62, y, 'VCSEL', ha='center', va='center', fontsize=5.8)
    axa.plot([1.20, 1.62], [y, y], color='0.35', lw=1.0)
    axa.add_patch(plt.Circle((1.80, y), 0.18, fc='white', ec='0.35', lw=0.8))
    axa.plot([1.98, 2.25], [y, y], color='0.35', lw=1.0)
    axa.plot([1.80, 1.80], [y - 0.18, y - 0.62], color='0.35', lw=0.8)
    axa.text(1.80, y - 0.78, 'PD', ha='center', va='center', fontsize=5.4,
             color='0.35')


def grating(x, y, col, h=0.20):
    for d in (-0.045, 0.0, 0.045):
        axa.plot([x + d, x + d], [y - h, y + h], color=col, lw=0.9)


# --- inline -----------------------------------------------------------
y1 = 3.95
source(y1)
axa.plot([2.25, 10.3], [y1, y1], color='0.35', lw=1.4)
axa.add_patch(Rectangle((2.45, y1 - 0.42), 1.55, 0.84, fc='#eaf1fb',
                        ec=BLUE, lw=0.7, ls=(0, (3, 2))))
for x in (2.75, 3.22, 3.69):
    grating(x, y1, BLUE)
axa.text(3.22, y1 + 0.60, 'references, stabilized', ha='center',
         fontsize=5.6, color=BLUE)
for x, c in ((5.1, VERM), (6.3, ORAN), (7.5, GREE), (8.7, PURP)):
    grating(x, y1, c)
axa.text(6.9, y1 - 0.62, 'sensors', ha='center', fontsize=5.6, color='0.35')
axa.text(0.05, y1 + 0.92, 'inline: every sensor looks through all three',
         fontsize=6.2, color='0.2')
for x in (5.1, 6.3, 7.5, 8.7):
    axa.add_patch(FancyArrowPatch((4.05, y1 + 0.30), (x, y1 + 0.30),
                                  arrowstyle='-', lw=0.4, color=VERM,
                                  alpha=0.55,
                                  connectionstyle='arc3,rad=-0.30'))

# --- branch -----------------------------------------------------------
y2 = 1.72
source(y2)
axa.plot([2.25, 10.3], [y2, y2], color='0.35', lw=1.4)
axa.plot([3.05, 3.05], [y2, y2 - 0.72], color='0.35', lw=1.0)
axa.plot([3.05, 4.60], [y2 - 0.72, y2 - 0.72], color='0.35', lw=1.4)
axa.add_patch(Rectangle((2.98, y2 - 0.10), 0.16, 0.20, fc='0.35', ec='none'))
axa.text(2.72, y2 + 0.26, 'coupler', fontsize=5.4, color='0.35', ha='center')
axa.add_patch(Rectangle((3.20, y2 - 1.04), 1.30, 0.64, fc='#eaf1fb',
                        ec=BLUE, lw=0.7, ls=(0, (3, 2))))
for x in (3.50, 3.85, 4.20):
    grating(x, y2 - 0.72, BLUE, h=0.16)
axa.text(3.85, y2 - 1.22, 'references, own stub', ha='center', fontsize=5.6,
         color=BLUE)
for x, c in ((5.1, VERM), (6.3, ORAN), (7.5, GREE), (8.7, PURP)):
    grating(x, y2, c)
axa.text(6.9, y2 + 0.42, 'sensors, nothing in front of them', ha='center',
         fontsize=5.6, color='0.35')
axa.text(0.05, y2 + 0.92, 'branch: the reference path is separate',
         fontsize=6.2, color='0.2')

FS.letter(axa, 'a')

# ==========================================================================
# (b) what the three references report when they sit inline
# ==========================================================================
grid = np.linspace(-210, 210, 300)
axb.plot(grid, axis_error(grid), color='0.35', lw=1.2,
         label='true axis error')
read10 = np.array([axis_error(r) + sum(law_a(REFS[j] - r, 0.10)
                                       for j in range(i))
                   for i, r in enumerate(REFS)])
fit = np.polyfit(REFS, read10, 2)
axb.fill_between(grid, axis_error(grid), np.polyval(fit, grid),
                 color=VERM, alpha=0.13, lw=0,
                 label='gap: added to every sensor')
axb.plot(grid, np.polyval(fit, grid), color=VERM, lw=1.0, ls=(0, (3, 2)),
         label='fit through the readings')
axb.plot(REFS, read10, 'o', color=VERM, ms=4.2,
         label='what each reference reads')
for r, v in zip(REFS, read10):
    if abs(v - axis_error(r)) > 0.3:
        axb.annotate('', xy=(r, v), xytext=(r, axis_error(r)),
                     arrowprops=dict(arrowstyle='-|>', color=VERM, lw=0.7,
                                     mutation_scale=6, shrinkA=0,
                                     shrinkB=1.5))
axb.text(-172, -6.6, '1st: nothing in front,\nreads the true error',
         fontsize=5.2, color='0.25', ha='left', va='top')
axb.text(14, 6.6, '2nd: shadowed by the 1st,\nreads $\\delta\\lambda_{k\\leftarrow j}=7.5$ pm too high',
         fontsize=5.2, color='0.25', ha='left', va='top')
axb.text(168, 11.6, '3rd: shadowed by two,\nreads 8.4 pm too high',
         fontsize=5.2, color='0.25', ha='right', va='top')
axb.set_xlim(-215, 215)
axb.set_ylim(-9, 17.5)
axb.set_xticks([-180, 0, 180])
axb.set_xlabel('band position [pm]', labelpad=1.5)
axb.set_ylabel('reported axis error [pm]', labelpad=1.5)
FS.letter(axb, 'b')
axb.legend(fontsize=5.2, loc='upper left', frameon=True, handlelength=1.5,
           labelspacing=0.18, borderaxespad=0.2)

# ==========================================================================
# (c) what survives the correction
# ==========================================================================
rr = np.array([0.003, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10])
inline_rms = [np.sqrt((residual(r) ** 2).mean()) for r in rr]
axc.semilogx(rr * 100, inline_rms, 'o-', color=VERM, ms=3.0,
             label='inline references')
axc.axhline(np.sqrt((axis_error(SENS) ** 2).mean()), color='0.35',
            ls=(0, (4, 2)), lw=0.9, label='no correction at all')
axc.axhline(0.02, color=BLUE, lw=1.4, label='references on a stub')
axc.set_xlim(0.25, 13)
axc.set_ylim(0, 7.5)
axc.set_xlabel('reference reflectivity [%]', labelpad=1.5)
axc.set_ylabel('RMS sensor error [pm]', labelpad=1.5)
FS.letter(axc, 'c')
axc.legend(fontsize=5.2, loc='upper left', frameon=True, handlelength=1.6,
           labelspacing=0.2, borderaxespad=0.25)
axc.grid(False, which='both', alpha=0.22)

os.makedirs('figs', exist_ok=True)
fig.savefig('figs/fig_s38_references.pdf', bbox_inches='tight',
            pad_inches=0.01)
fig.savefig('figs/fig_s38_references.png', dpi=200, bbox_inches='tight',
            pad_inches=0.01)
plt.close(fig)

print('reference readings at R = 10%%: %s pm'
      % ', '.join('%+.2f' % v for v in read10))
print('true axis error there:        %s pm'
      % ', '.join('%+.2f' % v for v in axis_error(REFS)))
print()
print('%-34s %8s %8s' % ('architecture', 'peak', 'RMS'))
print('%-34s %8.2f %8.2f' % ('no correction',
                             np.abs(axis_error(SENS)).max(),
                             np.sqrt((axis_error(SENS) ** 2).mean())))
for r in (0.10, 0.03, 0.01):
    res = residual(r)
    print('%-34s %8.2f %8.2f' % ('inline, R = %.0f%%' % (r * 100),
                                 np.abs(res).max(),
                                 np.sqrt((res ** 2).mean())))
print('%-34s %8.2f %8.2f' % ('on a stub, any R', 0.0, 0.0))
print()
print('saved figs/fig_s38_references.pdf')
