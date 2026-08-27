"""s29_ghostpath.py - the third-order ghost path, v2.

One CONTINUOUS light path with rounded U-turns at every reflection, its
width and opacity dropping by R at each bounce, next to the strong direct
echo of the same grating. Below, the time axis: the three direct echoes and
the weak ghost pulse arriving at tau_a - tau_b + tau_c, deeper than any of
its parents - a phantom grating, which is why it ghosts.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import figstyle as FS

FS.apply()

VERM, ORAN, GREE, BLUE = '#C74E0A', '#E8A200', '#009E73', '#0072B2'
GHOST = VERM

XB, XC, XA = 2.6, 5.4, 7.6            # b upstream; ghost of (a,b,c)
XIN = 0.5
FIBY = 4.0

fig, ax = plt.subplots(2, 1, figsize=(3.45, 2.35),
                       gridspec_kw=dict(height_ratios=[1.35, 0.85]))

# ---------------- geometry: one continuous path ----------------------------
g = ax[0]
g.plot([XIN - 0.2, 9.6], [FIBY, FIBY], color='0.62', lw=1.8,
       solid_capstyle='round')
for x0, col, lab in [(XB, ORAN, 'b'), (XC, GREE, 'c'), (XA, VERM, 'a')]:
    for dx in (-0.06, 0.0, 0.06):
        g.plot([x0 + dx, x0 + dx], [FIBY - 0.34, FIBY + 0.34], color=col,
               lw=1.6)
    g.text(x0, FIBY + 0.55, '$%s$' % lab, fontsize=8.5, ha='center',
           color=col)


def uturn(x, y1, y2, side, col, lw, al):
    """Rounded half-circle turn at x from lane y1 to lane y2."""
    r = abs(y2 - y1) / 2
    th = np.linspace(-np.pi / 2, np.pi / 2, 40)
    xc = x + (0.28 if side == 'right' else -0.28) * 0
    xs = x + (r * np.cos(th)) * (1 if side == 'right' else -1) * 0.55
    ys = (y1 + y2) / 2 + r * np.sin(th)
    g.plot(xs, ys, color=col, lw=lw, alpha=al, solid_capstyle='round')


def seg(x0, x1, y, col, lw, al, arrow=True):
    g.plot([x0, x1], [y, y], color=col, lw=lw, alpha=al,
           solid_capstyle='round')
    if arrow:
        xm = (x0 + x1) / 2
        d = 0.14 if x1 > x0 else -0.14
        g.annotate('', xy=(xm + d, y), xytext=(xm - d, y),
                   arrowprops=dict(arrowstyle='-|>', lw=0, color=col,
                                   alpha=al, mutation_scale=9))


# ghost path lanes (below fiber), width/alpha drop at each reflection
L1, L2, L3, L4 = 3.15, 2.45, 1.75, 1.05
seg(XIN, XA - 0.10, L1, GHOST, 2.4, 0.95)          # in, to a
uturn(XA - 0.05, L1, L2, 'right', GHOST, 1.6, 0.75)
seg(XA - 0.10, XB + 0.10, L2, GHOST, 1.6, 0.75)    # back to b
uturn(XB + 0.05, L2, L3, 'left', GHOST, 1.1, 0.55)
seg(XB + 0.10, XC - 0.10, L3, GHOST, 1.1, 0.55)    # forward to c
uturn(XC - 0.05, L3, L4, 'right', GHOST, 0.8, 0.40)
seg(XC - 0.10, XIN, L4, GHOST, 0.8, 0.40)          # out
for x0, y in [(XA, (L1 + L2) / 2), (XB, (L2 + L3) / 2),
              (XC, (L3 + L4) / 2)]:
    g.plot([x0, x0], [y - 0.02, FIBY - 0.36], color='0.85', lw=0.5,
           ls=(0, (1.5, 1.5)))
g.text(XIN - 0.12, L1, 'in', fontsize=6.4, ha='right', va='center',
       color='0.3')
g.text(XIN - 0.12, L4, 'ghost\nout', fontsize=6.0, ha='right', va='center',
       color=GHOST)
g.text(8.6, (L1 + L2) / 2, r'$\times R_a$', fontsize=6.4, color=VERM,
       va='center')
g.text(1.55, (L2 + L3) / 2, r'$\times R_b$', fontsize=6.4, color=ORAN,
       va='center')
g.text(6.35, (L3 + L4) / 2 + 0.02, r'$\times R_c$', fontsize=6.4,
       color=GREE, va='center')
g.set_xlim(-0.7, 10.2); g.set_ylim(0.55, 5.0)
g.axis('off')

# ---------------- time axis: why it ghosts ---------------------------------
t = ax[1]
TAUB, TAUC, TAUA = XB, XC, XA
TAUG = TAUA - TAUB + TAUC
for tau, col, h, lab in [(TAUB, ORAN, 1.0, r'$\tau_b$'),
                         (TAUC, GREE, 1.0, r'$\tau_c$'),
                         (TAUA, VERM, 1.0, r'$\tau_a$')]:
    t.plot([tau, tau], [0, h], color=col, lw=2.6, solid_capstyle='butt')
    t.text(tau, -0.34, lab, fontsize=7.5, ha='center', color=col)
t.plot([TAUG, TAUG], [0, 0.16], color=GHOST, lw=2.6,
       solid_capstyle='butt', alpha=0.75)
t.text(TAUG, -0.34, r'$\tau_g$', fontsize=7.5, ha='center', color=GHOST)
t.annotate('a phantom grating,\ndeeper than its parents',
           xy=(TAUG, 0.20), xytext=(9.15, 0.92), fontsize=6.2,
           ha='center', color=GHOST,
           arrowprops=dict(arrowstyle='-', lw=0.5, color=GHOST,
                           alpha=0.6))
t.text(0.35, 1.02, r'direct echoes $\propto R$', fontsize=6.4,
       color='0.25', va='top')
t.text(0.35, 0.55, r'ghost $\propto R^3$', fontsize=6.4, color=GHOST,
       va='top')
t.annotate('', xy=(10.9, 0), xytext=(0.2, 0),
           arrowprops=dict(arrowstyle='-|>', lw=0.7, color='0.45'))
t.text(10.85, 0.16, r'delay $\tau$', fontsize=6.6, ha='right',
       color='0.35')
t.set_xlim(-0.7, 11.2); t.set_ylim(-0.62, 1.25)
t.axis('off')

plt.tight_layout(h_pad=0.4)
fig.savefig('figs/fig_s29_ghostpath.png', dpi=140,
            bbox_inches='tight', pad_inches=0.03)
print('tau_g = %.1f (a=%.1f b=%.1f c=%.1f)' % (TAUG, TAUA, TAUB, TAUC))
print('saved figs/fig_s29_ghostpath.png')
