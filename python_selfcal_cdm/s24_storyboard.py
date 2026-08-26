"""s24_storyboard.py - the life of one code period, told as a picture.

One figure, one idea: the reader follows the code from the laser to the Bragg
wavelength. The trick that makes it legible is the shared horizontal scale. The
fiber is drawn on top with a DISTANCE axis; every time lane below uses the
matching TIME axis, time = 2nz/c, so the echo of each grating starts directly
underneath that grating, and the correlation peak lands there too. The dashed
verticals are therefore not decoration; they are the statement of the method.

All waveforms are computed with the real m-sequence, real delays and the real
correlator, not sketched. The lanes show the causal turn-on (echoes are silent
before their round trip); the correlation is evaluated in the periodic steady
state, one full period later.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle
import warnings; warnings.filterwarnings('ignore')
import common as C

plt.rcParams.update({'font.size': 8.5})

C_LIGHT, N_GROUP = 2.99792458e8, 1.468
B = 25e6
NBITS, N = 7, 127
OSR = 16
MSEQ = 1.0 - 2.0 * C._mls01(NBITS)
MSEQ = np.roll(MSEQ, -int(np.argmax(MSEQ == 1.0)))   # start on an "on" chip
DRIVE01 = 0.5 * (1.0 + MSEQ)              # on-off keying aligned with the code

z = np.array([16.3, 44.9, 73.4])          # metres
refl = np.array([0.10, 0.07, 0.05])
tau_chip = 2.0 * N_GROUP * z / C_LIGHT * B          # 4.0, 11.0, 18.0 chips
CH_SHOW = 26.0

# three periods so the causal start and the periodic steady state both exist
drive3 = np.tile(np.repeat(DRIVE01, OSR), 3)
t3 = np.arange(len(drive3)) / OSR                    # chips
rng = np.random.default_rng(5)
echo3 = []
for r, tc in zip(refl, tau_chip):
    sh = int(round(tc * OSR))
    e = np.zeros_like(drive3)
    e[sh:] = r * drive3[:len(drive3) - sh]           # causal, no wrap
    echo3.append(e)
rx3 = np.sum(echo3, axis=0) + rng.normal(0, 0.006, len(drive3))

# steady-state period for the correlator
rx_ss = rx3[N * OSR:2 * N * OSR].reshape(-1, OSR).mean(axis=1)
corr = C.periodic_xcorr(rx_ss - rx_ss.mean(), MSEQ)
lag = np.arange(N)

# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10.6, 6.6))
ax.set_xlim(-7.2, CH_SHOW + 0.6)
ax.set_ylim(-1.05, 11.9)
ax.axis('off')


def step(y, n, lines):
    ax.add_patch(Circle((-6.85, y), 0.40, fc='#1f4e79', ec='none', zorder=5))
    ax.text(-6.85, y, str(n), color='w', fontsize=8.5, ha='center', va='center',
            zorder=6, weight='bold')
    ax.text(-6.25, y, lines, fontsize=7.3, va='center', color='0.15')


sel3 = t3 < CH_SHOW

# --- fiber on top ------------------------------------------------------------
Y_FIB = 10.6
ax.plot([-1.2, CH_SHOW], [Y_FIB, Y_FIB], color='#555', lw=1.8)
ax.add_patch(Circle((-2.1, Y_FIB), 0.5, fc='#fdf3e7', ec='#1f4e79', lw=1.1))
ax.text(-2.1, Y_FIB, 'circ', fontsize=6.6, ha='center', va='center')
ax.add_patch(FancyBboxPatch((-7.1, Y_FIB - 0.55), 3.6, 1.1,
                            boxstyle='round,pad=0.05,rounding_size=0.1',
                            fc='#eaf1fb', ec='#1f4e79', lw=1.1))
ax.text(-5.3, Y_FIB, 'swept VCSEL\n+ code driver', fontsize=7.2, ha='center',
        va='center')
for zc, tc, nm in zip(z, tau_chip, ['FBG 1', 'FBG 2', 'FBG 3']):
    ax.add_patch(Rectangle((tc - 0.28, Y_FIB - 0.30), 0.56, 0.60, fc='#c0392b',
                           ec='#7b241c', hatch='///', lw=0.9, zorder=4))
    ax.text(tc, Y_FIB + 0.55, '%s\n%.0f m' % (nm, zc), fontsize=7, ha='center')
    ax.plot([tc, tc], [-0.25, Y_FIB - 0.32], color='0.75', ls='--', lw=0.7,
            zorder=0)
ax.annotate('', xy=(CH_SHOW - 0.2, Y_FIB - 0.80), xytext=(0.0, Y_FIB - 0.80),
            arrowprops=dict(arrowstyle='-|>', lw=0.9, color='#555'))
ax.text(CH_SHOW / 2, Y_FIB - 1.13,
        r'distance along the fiber, drawn to the time scale below:  $t = 2nz/c$',
        fontsize=7, ha='center', color='#555')

# --- lane 1: transmitted code ------------------------------------------------
Y1 = 8.35
step(Y1 + 0.30, 1, 'the sequence c(t)\non the laser current')
ax.plot(t3[sel3], Y1 + 0.62 * drive3[sel3], color='#1f4e79', lw=1.0)
ax.annotate('', xy=(2.0, Y1 + 0.88), xytext=(1.0, Y1 + 0.88),
            arrowprops=dict(arrowstyle='<->', lw=0.7))
ax.text(1.5, Y1 + 1.06, r'$T_{\rm chip}=1/B$', fontsize=6.4, ha='center')

# --- lanes 2: echoes ---------------------------------------------------------
Y2 = 6.95
step(Y2 - 0.55, 2, 'each grating returns the\nsame sequence, delayed by\n'
                   r'$\tau_k=2nz_k/c$, scaled by $R_k(\lambda)$')
cols = ['#c0392b', '#e67e22', '#b7950b']
for i, (e, tc, col) in enumerate(zip(echo3, tau_chip, cols)):
    y0 = Y2 - i * 0.80
    ax.plot(t3[sel3], y0 + 4.4 * e[sel3], color=col, lw=0.9)
    ax.plot(tc, y0 - 0.06, marker='^', ms=4, color=col, clip_on=False)

# --- lane 3: detector sum ----------------------------------------------------
Y3 = 3.95
step(Y3 + 0.25, 3, 'one photodiode sees\nthe sum, plus noise')
ax.plot(t3[sel3], Y3 + 2.6 * rx3[sel3], color='#7d3c98', lw=0.9)

# --- correlator block --------------------------------------------------------
Y4 = 2.42
ax.add_patch(FancyBboxPatch((7.6, Y4 - 0.34), 10.6, 0.72,
                            boxstyle='round,pad=0.05,rounding_size=0.1',
                            fc='#eafaf0', ec='#1f4e79', lw=1.1))
ax.text(12.9, Y4, 'correlate with the transmitted sequence', fontsize=7.4,
        ha='center', va='center')
step(Y4, 4, 'matched filter')
ax.annotate('', xy=(12.9, Y4 + 0.42), xytext=(12.9, Y3 - 0.45),
            arrowprops=dict(arrowstyle='-|>', lw=1.0, color='#333'))
ax.annotate('', xy=(12.9, 1.72), xytext=(12.9, Y4 - 0.40),
            arrowprops=dict(arrowstyle='-|>', lw=1.0, color='#333'))

# --- lane 4: correlation output ---------------------------------------------
Y5 = 0.42
step(Y5 + 0.62, 5, 'one peak per grating:\nlag = position, height = \n'
                   r'$R_k$ at the current sweep $\lambda$')
m = lag <= CH_SHOW
ax.plot(lag[m], Y5 + np.maximum(corr[m], -0.004) / corr.max() * 1.1,
        color='#2980b9', lw=1.2)
for tc, col in zip(tau_chip, cols):
    k = int(round(tc))
    ax.plot(k, Y5 + corr[k] / corr.max() * 1.1 + 0.10, 'v', ms=5, color=col)
ax.annotate('', xy=(CH_SHOW - 0.2, Y5 - 0.42), xytext=(0.0, Y5 - 0.42),
            arrowprops=dict(arrowstyle='-|>', lw=0.9, color='#555'))
ax.text(CH_SHOW / 2, Y5 - 0.72, 'lag [chips]  =  round-trip delay  =  position',
        fontsize=7, ha='center', color='#555')

# --- step 6: the sweep builds the map ---------------------------------------
ax.add_patch(FancyBboxPatch((18.6, 1.35), 7.6, 1.75,
                            boxstyle='round,pad=0.12,rounding_size=0.12',
                            fc='#eaf1fb', ec='#1f4e79', lw=0.9))
ax.add_patch(Circle((19.3, 2.72), 0.36, fc='#1f4e79', ec='none', zorder=5))
ax.text(19.3, 2.72, '6', color='w', fontsize=8, ha='center', va='center',
        zorder=6, weight='bold')
ax.text(19.95, 2.20, 'repeat at every sweep wavelength:\nthe peak heights trace'
                     ' each grating\'s\nspectrum (the map of Fig. 1b)',
        fontsize=6.8, va='center', color='#1f4e79')

fig.tight_layout()
fig.savefig('figs/fig_s24_storyboard.png', dpi=150, bbox_inches='tight')
print('gratings at %s m -> delays %s chips'
      % (z.tolist(), np.round(tau_chip, 1).tolist()))
print('correlation peaks at the grating bins: %s (heights %s)'
      % ([int(round(tc)) for tc in tau_chip],
         np.round([corr[int(round(tc))] for tc in tau_chip], 4).tolist()))
print('saved figs/fig_s24_storyboard.png')
