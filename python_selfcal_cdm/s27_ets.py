"""s27_ets.py - equivalent-time sampling in one small picture.

Top: the periodic coded probe over several sequence periods; a slow
converter takes one sparse comb of samples per period, phase-stepped by the
programmable delay. Bottom: the samples reassembled into one dense record
of a single period. Single-column figure for Sec. Acquisition.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS

FS.apply()

NSH = 24                       # chips shown per period
NPER = 3                       # periods shown
NS = 4                         # samples per chip after reassembly
_MS = 1.0 - 2.0 * C._mls01(7)
code = 0.5 * (1 + _MS[:NSH])
OS = 24
sig1 = np.repeat(code, OS)
t1 = np.arange(len(sig1)) / OS
sig = np.tile(sig1, NPER)
t = np.arange(len(sig)) / OS

fig, ax = plt.subplots(2, 1, figsize=(3.45, 1.9),
                       gridspec_kw=dict(height_ratios=[1.15, 0.85]))

a = ax[0]
a.plot(t, sig, color='0.45', lw=0.7)
COLS = [FS.BLUE, FS.GREEN, FS.VERM]
NSAMP = 6                      # sparse samples per period
for p in range(NPER):
    ph = p / (NPER * 1.0)      # phase step between periods
    ts = (np.arange(NSAMP) + 0.12 + ph) * (NSH / NSAMP) + p * NSH
    ts = ts[ts < (p + 1) * NSH]
    a.plot(ts, np.interp(ts, t, sig), ls='none', marker='o', ms=3.2,
           color=COLS[p], mec='white', mew=0.4, zorder=5)
    a.axvline(p * NSH, color='0.8', lw=0.6)
    a.text(p * NSH + NSH / 2, 1.22, 'period %d' % (p + 1), fontsize=6.4,
           ha='center', color=COLS[p])
a.annotate('', xy=(NSH * 0.62, -0.28), xytext=(NSH * 0.12, -0.28),
           annotation_clip=False,
           arrowprops=dict(arrowstyle='-|>', lw=0.6, color='0.4'))
a.text(NSH * 0.37, -0.44, 'delay step', fontsize=6.2, ha='center',
       color='0.4')
a.set_xlim(0, NPER * NSH); a.set_ylim(-0.15, 1.35)
a.set_xticks([]); a.set_yticks([])
for sp in a.spines.values():
    sp.set_visible(False)
a.set_title('slow converter, one sparse comb per period', fontsize=7.5)

b = ax[1]
b.plot(t1, sig1, color='0.45', lw=0.7)
td = np.arange(0.12, NSH, 1.0 / NS)
allc = np.concatenate([[COLS[p]] * 1 for p in range(NPER)])
for i, tsamp in enumerate(td):
    b.plot([tsamp], [np.interp(tsamp, t1, sig1)], ls='none', marker='o',
           ms=2.0, color=COLS[i % NPER], mec='none', zorder=5)
b.set_xlim(0, NSH); b.set_ylim(-0.15, 1.3)
b.set_xticks([]); b.set_yticks([])
for sp in b.spines.values():
    sp.set_visible(False)
b.set_title('reassembled: one dense record of one period', fontsize=7.5)
b.set_xlabel(r'record time $N n_s/f_s$ (147 $\mu$s at 3.6 MS/s)',
             fontsize=7)

plt.tight_layout()
fig.savefig('figs/fig_s27_ets.png', dpi=140)
print('saved figs/fig_s27_ets.png')
