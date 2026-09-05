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

The branch case assumes ideal reference readings, isolated from the sensors.
Mutual shadowing among the serial reference gratings is not modeled in that
case. The zero residual is therefore an assumption of the comparison, not a
prediction for an arbitrary serial reference branch. The extra 6 dB of optical
loss reduces the detector-noise SNR from 99 to 87 dB in the 20log10 convention.

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
axa.set_xlim(0, 10.8)
axa.set_ylim(0, 7.4)
axa.set_aspect('equal', adjustable='box')
axa.axis('off')

INK = '#40464D'
SENSOR = VERM


def fiber(xs, ys):
    axa.plot(xs, ys, color=INK, lw=1.05, solid_capstyle='round', zorder=1)


def component(x, y, width, height, label):
    axa.add_patch(Rectangle((x, y - height / 2), width, height,
                            facecolor='white', edgecolor=INK, lw=0.8, zorder=3))
    axa.text(x + width / 2, y, label, ha='center', va='center',
             fontsize=6.2, color=INK, zorder=4)


def source(y):
    component(0.08, y, 1.28, 0.62, 'VCSEL')
    fiber([1.36, 1.88], [y, y])
    fiber([2.42, 2.95], [y, y])
    fiber([2.15, 2.15], [y - 0.27, y - 0.80])
    axa.add_patch(plt.Circle((2.15, y), 0.27, facecolor='white',
                            edgecolor=INK, lw=0.8, zorder=3))
    axa.add_patch(FancyArrowPatch((2.00, y + 0.10), (2.29, y - 0.10),
                                  connectionstyle='arc3,rad=-0.65',
                                  arrowstyle='-|>', mutation_scale=4,
                                  lw=0.65, color=INK, zorder=4))
    axa.text(2.15, y + 0.55, 'Circulator', ha='center', fontsize=5.6,
             color=INK)
    component(1.77, y - 1.03, 0.76, 0.46, 'PD')


def grating(x, y, color):
    axa.add_patch(Rectangle((x - 0.11, y - 0.24), 0.22, 0.48,
                            facecolor='white', edgecolor=color,
                            lw=0.75, zorder=3))
    for offset in (-0.15, 0.0, 0.15):
        axa.plot([x - 0.09, x + 0.09],
                 [y + offset - 0.055, y + offset + 0.055],
                 color=color, lw=0.65, zorder=4)


def references(xs, y, label_above=True):
    left, right = xs[0] - 0.36, xs[-1] + 0.36
    axa.add_patch(Rectangle((left, y - 0.41), right - left, 0.82,
                            facecolor='#EDF5FA', edgecolor=BLUE, lw=0.65,
                            linestyle=(0, (3, 2)), zorder=0))
    for i, x in enumerate(xs, 1):
        grating(x, y, BLUE)
        axa.text(x, y - 0.69, '$R_%d$' % i, fontsize=5.8,
                 color=BLUE, ha='center', va='center')
    label_y = y + 0.60 if label_above else y - 1.12
    axa.text((left + right) / 2, label_y, 'References',
             ha='center', fontsize=6.0, color=BLUE)


def sensors(y):
    for i, x in enumerate((6.25, 7.35, 8.45, 9.55), 1):
        grating(x, y, SENSOR)
        axa.text(x, y - 0.69, '$S_%d$' % i, fontsize=5.8,
                 color=INK, ha='center', va='center')
    axa.text(7.90, y + 0.60, 'Sensors', ha='center', fontsize=6.0,
             color=INK)


# Keep corresponding components aligned between the two arrangements.
upper, lower = 5.95, 2.65
axa.text(0.08, 7.10, 'Inline references', fontsize=7.0,
         fontweight='bold', color=INK)
source(upper)
fiber([2.95, 10.40], [upper, upper])
references((3.65, 4.30, 4.95), upper)
sensors(upper)

axa.plot([0.08, 10.40], [4.30, 4.30], color='#DDE1E5', lw=0.55)
axa.text(0.08, 3.83, 'Separate reference branch', fontsize=7.0,
         fontweight='bold', color=INK)
source(lower)
fiber([2.95, 10.40], [lower, lower])
sensors(lower)

# A boxed splitter makes the extra optical component explicit.
split_x, ref_y = 3.48, 1.22
axa.add_patch(Rectangle((3.06, lower - 0.30), 0.84, 0.60,
                        facecolor='white', edgecolor=INK, lw=0.8, zorder=3))
axa.plot([3.06, 3.90], [lower, lower], color=INK, lw=0.7, zorder=4)
axa.plot([split_x, split_x], [lower, lower - 0.30],
         color=INK, lw=0.7, zorder=4)
fiber([split_x, split_x, 5.55], [lower - 0.30, ref_y, ref_y])
axa.text(split_x, lower + 0.55, 'Coupler', ha='center',
         fontsize=5.6, color=INK)
references((3.90, 4.50, 5.10), ref_y, label_above=False)

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
axb.text(-172, -6.6, '$R_1$: no bias',
         fontsize=5.2, color='0.25', ha='left', va='top')
axb.text(14, 6.6, '$R_2$: +7.5 pm bias',
         fontsize=5.2, color='0.25', ha='left', va='top')
axb.text(168, 11.6, '$R_3$: +8.4 pm bias',
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
