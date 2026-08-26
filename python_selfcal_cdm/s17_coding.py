"""s17_coding.py - what the code actually does, in the time domain.

The explanatory figure that the earlier draft was missing. Three panels, all
computed with the real sequence rather than sketched:

  (a) the drive: one period of a maximal-length sequence, on-off keying the
      VCSEL current;
  (b) what the detector sees: the same sequence returned three times, once per
      grating, each delayed by its own round trip and scaled by its reflectance,
      summed and buried in noise. Nothing in that trace identifies a grating;
  (c) what the correlator recovers: correlating with the transmitted sequence
      collapses each return into a peak at its own delay, on a floor set by the
      periodic autocorrelation side lobe. The delay axis is the position axis,
      and the peak height at that lag is the reflectance at the wavelength the
      sweep happened to be at.

The inset in (c) shows the shape of a single peak: the triangular main lobe of a
rectangular-chip sequence, two chips wide, which is what ties the spatial
resolution to the chip rate.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import common as C

plt.rcParams.update({'font.size': 8.5})

NBITS = 7
N = 2 ** NBITS - 1
OSR = 8                       # samples per chip in the time-domain plot
B = 25e6                      # chip rate [chip/s]
C_LIGHT, N_GROUP = 2.99792458e8, 1.468

seq = 1.0 - 2.0 * C._mls01(NBITS)       # bipolar sequence
drive01 = 0.5 * (1.0 + seq)             # on-off keyed laser current
drive = np.repeat(drive01, OSR).astype(float)
t_us = np.arange(len(drive)) / (B * OSR) * 1e6

# three gratings, at 4.1, 9.6 and 15.2 m, with different reflectances
z = np.array([41.0, 102.0, 184.0])
tau_chip = 2.0 * N_GROUP * z / C_LIGHT * B
refl = np.array([0.10, 0.07, 0.05])

rng = np.random.default_rng(3)
rx = np.zeros_like(drive)
for zi, ri in zip(tau_chip, refl):
    shift = int(round(zi * OSR))
    rx += ri * np.roll(drive, shift)
rx = rx + rng.normal(0, 0.006, rx.shape)

# correlate against the bipolar sequence, one sample per chip
rx_chip = rx.reshape(-1, OSR).mean(axis=1)
corr = C.periodic_xcorr(rx_chip - rx_chip.mean(), seq)
lag = np.arange(N)

fig, ax = plt.subplots(1, 3, figsize=(10.6, 3.2))

# --- (a) the drive ---------------------------------------------------------
n_show = 24 * OSR
ax[0].step(t_us[:n_show], drive[:n_show], where='post', color='#1f4e79', lw=1.2)
ax[0].set_ylim(-0.25, 1.35)
ax[0].set_xlabel(r'time [$\mu$s]'); ax[0].set_ylabel('drive (on / off)')
ax[0].set_title('(a) The code on the laser current', fontsize=9)
ax[0].text(0.03, 0.30, 'm-sequence, N = %d, %.0f Mchip/s' % (N, B / 1e6),
           transform=ax[0].transAxes, fontsize=7, color='#1f4e79')
i0 = int(np.argmax(drive[:n_show] > 0.5))
ax[0].annotate('', xy=(t_us[i0], 1.16), xytext=(t_us[i0 + OSR], 1.16),
               arrowprops=dict(arrowstyle='<->', lw=0.8))
ax[0].text(t_us[i0 + OSR // 2], 1.20, r'$T_{\rm chip}$', ha='center', fontsize=7)
ax[0].grid(True, alpha=0.2)

# --- (b) what the detector sees --------------------------------------------
ax[1].plot(t_us[:n_show], rx[:n_show], color='#7d3c98', lw=1.0)
ax[1].set_xlabel(r'time [$\mu$s]'); ax[1].set_ylabel('detected power [a.u.]')
ax[1].set_title('(b) Three gratings, summed and noisy', fontsize=9)
ax[1].text(0.02, 0.92, 'no grating is identifiable here',
           transform=ax[1].transAxes, fontsize=7, color='#7d3c98')
ax[1].grid(True, alpha=0.2)

# --- (c) after correlation --------------------------------------------------
ax[2].plot(lag, corr, color='#2980b9', lw=1.0)
for zi, ri, zz in zip(tau_chip, refl, z):
    ax[2].annotate('%.1f m' % zz, xy=(zi, np.interp(zi, lag, corr)),
                   xytext=(zi, np.interp(zi, lag, corr) + 0.011),
                   ha='center', fontsize=7, color='#c0392b',
                   arrowprops=dict(arrowstyle='-', lw=0.6, color='#c0392b'))
ax[2].set_xlabel('lag [chip]  ->  round-trip delay  ->  position')
ax[2].set_ylabel('correlator output')
ax[2].set_title('(c) After despreading', fontsize=9)
ax[2].set_xlim(0, 60)
ax[2].grid(True, alpha=0.2)

axin = ax[2].inset_axes([0.58, 0.14, 0.39, 0.40])
k0 = int(round(tau_chip[0]))
sel = slice(k0 - 4, k0 + 5)
axin.plot(lag[sel] - k0, corr[sel], 'o-', color='#2980b9', ms=3, lw=1.0)
axin.set_title('peak shape', fontsize=6.5, pad=2)
axin.tick_params(labelsize=6)
axin.set_xlabel('lag [chip]', fontsize=6, labelpad=1)
axin.text(0.05, 0.72, '2 chips wide', transform=axin.transAxes, fontsize=6,
          color='0.35')

fig.tight_layout()
fig.savefig('figs/fig_s17_coding.png', dpi=150, bbox_inches='tight')

print('gratings at %s m -> lags %s chip' % (z.tolist(), np.round(tau_chip, 1).tolist()))
print('recovered peak heights: %s (true reflectances %s)'
      % (np.round([corr[int(round(t))] for t in tau_chip], 4).tolist(), refl.tolist()))
mask = np.ones(N, bool)
for t in tau_chip:
    k = int(round(t)); mask[max(0, k - 2):k + 3] = False
print('floor away from the peaks: %.4f, largest peak %.4f, ratio %.3f'
      % (np.abs(corr[mask]).max(), corr.max(), np.abs(corr[mask]).max() / corr.max()))
print('note: the unipolar drive halves the recovered amplitude, as in the link budget')
print('saved figs/fig_s17_coding.png')
