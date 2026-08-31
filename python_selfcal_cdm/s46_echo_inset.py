"""s46_echo_inset.py - the echo inset of the measurement-chain figure.

Version of 2026-08-31, replacing the inset_echosum drawn by s24_concept v11
(commit 354aefb). Same code, same delays and same weights as before, so the
picture stays consistent with the drive inset (s35) and with the delay
positions of the correlated row (s37). What changed: every trace is named
(FBG1 .. FBG4 with its term A_k(lambda_m) c(t - tau_k), and the sum P_m(t)),
and the two axes are drawn and labelled, time along x and optical power
along y.

The code is the m-sequence of length 127 rolled to start with a "1" chip.
The drive inset shows its first ten chips at every wavelength step, the echo
inset its first 30 chip periods after the transmit instant t = 0. Each
grating returns the code delayed by tau_k and weighted by A_k(lambda_m), the
amplitude that reaches the detector. The photodiode sees only the sum.

Output: figs/fig1_panels/inset_echosum.pdf (and .png, .svg), consumed by
Draft_v2/fig1_concept.tex at a width of 4.8 cm. Drawn at that width, so the
font sizes below are the printed ones.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import common as C
import figstyle as FS

FS.apply()

ECOL = ['#C74E0A', '#E8A200', '#009E73', '#CC79A7']   # grating colours, as in
                                                      # the rest of the figure
FS_LBL = 5.8
FS_SML = 5.0

# --- the code, identical construction to s35 and to s24 v11 --------------
_MS = 1.0 - 2.0 * C._mls01(7)
_MS = np.roll(_MS, -int(np.argmax(_MS == 1.0)))
OS = 6
drv = np.repeat(0.5 * (1 + _MS[:34]), OS)
tt = np.arange(len(drv)) / float(OS)
TW = 30.0                                  # chip periods shown
W = tt < TW

taus = [2.2, 6.0, 10.2, 14.2]              # delays in chips
amps = [1.00, 0.74, 0.52, 0.36]            # A_k(lambda_m), relative
offs = (5.0, 4.0, 3.0, 2.0)                # trace baselines
H = 0.85                                   # trace height
SUM = 0.42                                 # sum scale

fig = plt.figure(figsize=(1.89, 1.20))     # 4.8 cm wide
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])

ech = []
for k, (tau, a, col, off) in enumerate(zip(taus, amps, ECOL, offs)):
    e = np.zeros_like(drv)
    sh = int(tau * OS)
    e[sh:] = a * drv[:len(drv) - sh]
    ech.append(e)
    ax.plot(tt[W], off + H * e[W], color=col, lw=0.7)
    ax.text(TW + 0.6, off + 0.5 * H * a, r'FBG%d: $A_%d(\lambda_m)\,c(t-\tau_%d)$'
            % (k + 1, k + 1, k + 1), fontsize=FS_SML, color=col, ha='left',
            va='center')
sm = np.sum(ech, axis=0)
ax.plot(tt[W], SUM * sm[W], color='0.1', lw=0.8)
ax.text(TW + 0.6, SUM * sm[W].max() * 0.5, r'sum: $P_m(t)$', fontsize=FS_SML,
        color='0.1', ha='left', va='center')

# transmit instant and the tau_2 dimension
ax.axvline(0, color='0.45', lw=0.5, ls=(0, (2, 1.5)))
ax.text(0.3, offs[0] + H + 0.32, r'$t=0$', fontsize=FS_SML, color='0.45',
        ha='left', va='bottom')
# tau_2 just under the FBG2 baseline, in the gap above the FBG3 trace
y2 = offs[1] - 0.28
ax.annotate('', xy=(taus[1], y2), xytext=(0, y2),
            arrowprops=dict(arrowstyle='-|>', lw=0.6, color=ECOL[1],
                            mutation_scale=5))
ax.text(taus[1] + 0.5, y2, r'$\tau_2$', fontsize=FS_LBL, color=ECOL[1],
        ha='left', va='center')

# axes: time along x, optical power along y
X0, Y0 = -1.6, -0.55
ax.annotate('', xy=(TW + 0.4, Y0), xytext=(X0, Y0),
            arrowprops=dict(arrowstyle='-|>', lw=0.7, color='0.3',
                            mutation_scale=6))
ax.text(TW + 0.6, Y0, r'$t$', fontsize=FS_LBL, color='0.3', ha='left',
        va='center')
ax.annotate('', xy=(X0, offs[0] + H + 0.75), xytext=(X0, Y0),
            arrowprops=dict(arrowstyle='-|>', lw=0.7, color='0.3',
                            mutation_scale=6))
ax.text(X0 - 0.3, offs[0] + H + 0.75, r'$P$ (optical power)', fontsize=FS_LBL,
        color='0.3', ha='left', va='bottom')

ax.set_xlim(X0 - 0.5, TW + 13.5)
ax.set_ylim(Y0 - 0.25, offs[0] + H + 1.35)
ax.axis('off')

os.makedirs('figs/fig1_panels', exist_ok=True)
fig.savefig('figs/fig1_panels/inset_echosum.pdf', bbox_inches='tight',
            pad_inches=0.01, transparent=True)
fig.savefig('figs/fig1_panels/inset_echosum.png', dpi=400, bbox_inches='tight',
            pad_inches=0.01)
fig.savefig('figs/fig1_panels/inset_echosum.svg', bbox_inches='tight',
            pad_inches=0.01, transparent=True)
plt.close(fig)
print('chips [:10]:', (0.5 * (1 + _MS))[:10].astype(int))
print('saved figs/fig1_panels/inset_echosum.pdf, .png, .svg')
