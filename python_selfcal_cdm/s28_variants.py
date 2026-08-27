"""s28_variants.py - the two coded-sweep variants in one small figure.

Top: synchronous variant (Markowski 2023): every wavelength carries its own
sequence, all M codes overlap in one observation window; cross-correlation
separates them, so the family must hold at least M sequences.
Bottom: stepped variant (Koo 1999, this paper): one sequence reused at every
wavelength step, one code on the fiber at a time; autocorrelation separates
the gratings, the frame is M times longer, the family size is irrelevant.
Single-column figure for Sec. II-B.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS

FS.apply()

COLS = ['#C74E0A', '#E8A200', '#009E73', '#CC79A7']
M = 4
NSH = 18
OS = 10
_MS = 1.0 - 2.0 * C._mls01(7)

fig, ax = plt.subplots(2, 1, figsize=(3.45, 2.05))

# --- synchronous: M codes at once ------------------------------------------
a = ax[0]
for m in range(M):
    code = 0.5 * (1 + _MS[7 * m:7 * m + NSH])       # different codes
    y = np.repeat(code, OS)
    t = np.arange(len(y)) / OS
    a.plot(t, (M - 1 - m) * 1.35 + y, color=COLS[m], lw=0.8)
    a.text(-0.5, (M - 1 - m) * 1.35 + 0.5, r'$\lambda_%d$' % (m + 1),
           fontsize=6.6, ha='right', va='center', color=COLS[m])
a.annotate('', xy=(NSH, -0.62), xytext=(0, -0.62), annotation_clip=False,
           arrowprops=dict(arrowstyle='<->', lw=0.6, color='0.4'))
a.text(NSH / 2, -1.25, 'one window: all $M$ codes overlap', fontsize=6.2,
       ha='center', color='0.35')
a.set_xlim(-1.6, NSH + 0.3); a.set_ylim(-0.8, M * 1.35 + 0.1)
a.axis('off')
a.set_title('synchronous: cross-correlation separates, family $\\geq M$',
            fontsize=7.2)

# --- stepped: one code per step --------------------------------------------
b = ax[1]
code = 0.5 * (1 + _MS[:NSH])
y1 = np.repeat(code, OS)
t1 = np.arange(len(y1)) / OS
for m in range(M):
    b.plot(t1 + m * NSH, 0.15 + y1, color=COLS[m], lw=0.8)
    b.text(m * NSH + NSH / 2, 1.55, r'$\lambda_%d$' % (m + 1), fontsize=6.6,
           ha='center', color=COLS[m])
    if m:
        b.axvline(m * NSH, color='0.85', lw=0.6)
b.annotate('', xy=(M * NSH, -0.62), xytext=(0, -0.62), annotation_clip=False,
           arrowprops=dict(arrowstyle='<->', lw=0.6, color='0.4'))
b.text(M * NSH / 2, -1.30, r'$M$ windows: the same code, reused',
       fontsize=6.2, ha='center', color='0.35')
b.set_xlim(-1.6 * M, M * NSH + 0.3); b.set_ylim(-0.85, 2.0)
b.axis('off')
b.set_title('stepped: autocorrelation separates, family $=1$',
            fontsize=7.2)

plt.tight_layout(h_pad=1.4)
fig.savefig('figs/fig_s28_variants.png', dpi=140,
            bbox_inches='tight', pad_inches=0.03)
print('saved figs/fig_s28_variants.png')
