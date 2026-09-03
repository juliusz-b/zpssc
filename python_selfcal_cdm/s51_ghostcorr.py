"""s51_ghostcorr.py - the delay profile of the three gratings of Fig. 6.

Time-domain companion of s50_ghostpath: three gratings read with a unipolar
m-sequence at 100 Mchip/s (one chip is 1.02 m of fibre). The photodiode record is built from every
path up to the third order (three direct returns of amplitude ~R, five
third-order ghosts of amplitude ~R^3, all with the (1-R) transmission of every
grating crossed on the way), the record mean is removed and the record is
correlated with the bipolar replica, as in the receiver of the paper. The
magnitude of the correlation is drawn against position on a logarithmic
scale, so that the ghosts, two orders of magnitude below the direct returns
at R = 10%, are visible next to the code side lobes at -1/N.

  (a) The spacing of Fig. 6 (z_b : z_c : z_a = 1 : 2.4 : 5): five ghosts on
      five distinct positions, one pair sharing the position z_a - z_b + z_c.
  (b) Uniform spacing: the ghost (c,b,c) lands exactly under grating a, and
      the three ghosts (a,b,c), (c,b,a), (a,c,a) pile up on one position.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import common as C
import figstyle as FS

FS.apply()

# ---------------------------------------------------------------------------
C_LIGHT, NG = 2.998e8, 1.468
CHIP_RATE = 100e6                 # chip rate of the designed column
SPC = 16                          # samples per chip
NBITS = 9                         # N = 511
R = 0.10
ZMAXPLOT = 44.0
LAYOUTS = {'a': {'b': 4.0, 'c': 9.6, 'a': 20.0},      # same ratios as Fig. 6
           'b': {'b': 4.0, 'c': 8.0, 'a': 12.0}}      # uniform 4 m
COL = {'a': FS.VERM, 'b': FS.ORANGE, 'c': FS.GREEN}
GREY = '0.45'


def tint(col, w):
    c = np.array(to_rgb(col)); g = np.array([0.78, 0.78, 0.78])
    return tuple(w * c + (1 - w) * g)


GH = {g: tint(COL[g], 0.45) for g in COL}
GHT = {g: tint(COL[g], 0.7) for g in COL}

M_PER_CHIP = C_LIGHT / NG / CHIP_RATE / 2.0   # round trip: one chip = this many metres
print('one chip = %.2f m of fibre' % M_PER_CHIP)


# ---------------------------------------------------------------------------
# every path up to the third order
# ---------------------------------------------------------------------------
def paths(Z):
    """[(sequence of gratings, amplitude, equivalent position), ...]"""
    order = sorted(Z, key=Z.get)

    def trans(z0, z1):
        t = 1.0
        for g in order:
            if min(z0, z1) < Z[g] < max(z0, z1):
                t *= (1.0 - R)
        return t

    def one(seq):
        amp, zz, length = 1.0, 0.0, 0.0
        for g in seq:
            amp *= trans(zz, Z[g]) * R
            length += abs(Z[g] - zz)
            zz = Z[g]
        amp *= trans(zz, 0.0)
        return amp, (length + zz) / 2.0

    out = [((g,),) + one((g,)) for g in order]
    for g1 in order:
        for g2 in order:
            for g3 in order:
                if Z[g2] < Z[g1] and Z[g2] < Z[g3]:
                    out.append(((g1, g2, g3),) + one((g1, g2, g3)))
    return out


def delay_profile(Z):
    code01 = C._mls01(NBITS)
    n = code01.size
    tx = np.repeat(code01.astype(float), SPC)
    rx = np.zeros_like(tx)
    for seq, amp, zpos in paths(Z):
        rx += amp * np.roll(tx, int(round(zpos / M_PER_CHIP * SPC)))
    replica = np.repeat(C._to_pm1(code01).astype(float), SPC)
    rx = rx - rx.mean()
    corr = np.fft.ifft(np.fft.fft(rx) * np.conj(np.fft.fft(replica))).real
    corr /= (n * SPC / 2.0)                      # a unit return gives a peak of 1
    corr -= np.median(corr)                      # constant offset of the periodic correlation, -1/N per return
    zaxis = np.arange(tx.size) / SPC * M_PER_CHIP
    return zaxis, corr, n


fig, axes = plt.subplots(2, 1, figsize=(3.45, 3.2), sharex=True,
                         gridspec_kw=dict(hspace=0.08, left=0.15, right=0.985, bottom=0.125, top=0.955))

for let, Z in LAYOUTS.items():
    ax = axes[0] if let == 'a' else axes[1]
    zaxis, corr, n = delay_profile(Z)
    P = paths(Z)
    ref = P[0][1]                                 # amplitude of the first direct return
    m = zaxis <= ZMAXPLOT
    ax.plot(zaxis[m], np.abs(corr[m]) / ref, color='0.25', lw=0.8, zorder=3)
    ax.axhline(1.0 / n, color='0.6', lw=0.7, ls='--', zorder=1)
    ax.text(ZMAXPLOT - 0.5, 1.0 / n * 1.25, r'$1/N$', color='0.45', ha='right', va='bottom', fontsize=6.4)
    # markers: one per distinct position, at the summed amplitude, ghosts in
    # the greyed colour of their first grating, direct returns in full colour
    pos = {}
    for seq, amp, zpos in P:
        pos.setdefault(round(zpos, 3), []).append((seq, amp))
    for zpos, lst in pos.items():
        tot = sum(a for _, a in lst) / ref
        direct = [s for s, _ in lst if len(s) == 1]
        ghosts = [s for s, _ in lst if len(s) == 3]
        col = COL[direct[0][0]] if direct else GH[ghosts[0][0]]
        ax.plot(zpos, tot, 'v', color=col, ms=3.8, mec='white', mew=0.4, zorder=5, clip_on=False)
        if direct:
            ax.text(zpos, 1.9, r'$%s$' % direct[0][0], color=COL[direct[0][0]], ha='center', va='bottom', fontsize=7)
    ax.set_yscale('log')
    ax.set_ylim(3e-4, 3.0)
    ax.set_yticks([1e-3, 1e-2, 1e-1, 1])
    FS.letter(ax, let)
    print('layout', let, {k: [s for s, _ in v] for k, v in pos.items()})

# ghost labels, panel (a): every ghost on its own position
Pa = paths(LAYOUTS['a']); ref = Pa[0][1]
LAB = {('c', 'b', 'c'): (r'$(c,b,c)$', 0), ('a', 'b', 'c'): (r'$(a,b,c)$ $(c,b,a)$', 0),
       ('a', 'c', 'a'): (r'$(a,c,a)$', 0), ('a', 'b', 'a'): (r'$(a,b,a)$', 0)}
for seq, amp, zpos in Pa:
    if seq in LAB:
        tot = amp / ref * (2.0 if seq == ('a', 'b', 'c') else 1.0)
        axes[0].text(zpos, tot * 2.3, LAB[seq][0], color=GHT[seq[0]], ha='center', va='bottom', fontsize=6.0)
axes[0].text(0.98, 0.94, 'spacing of Fig. 6', transform=axes[0].transAxes, ha='right', va='top', fontsize=6.6)

# ghost labels, panel (b): collisions
Pb = paths(LAYOUTS['b']); ref = Pb[0][1]
amp_cbc = [a for s, a, _ in Pb if s == ('c', 'b', 'c')][0]
axes[1].annotate(r'$(c,b,c)$ under $a$', xy=(12.0, 1.0), xytext=(17.5, 0.7),
                 color=GHT['c'], fontsize=6.0, ha='left', va='center',
                 arrowprops=dict(arrowstyle='-', color=GHT['c'], lw=0.6, shrinkA=0, shrinkB=3))
three = sum(a for s, a, _ in Pb if s in (('a', 'b', 'c'), ('c', 'b', 'a'), ('a', 'c', 'a'))) / ref
axes[1].text(17.3, three * 1.6, r'$(a,b,c)$ $(c,b,a)$ $(a,c,a)$', color=GHT['a'], ha='left', va='center', fontsize=6.0)
amp_aba = [a for s, a, _ in Pb if s == ('a', 'b', 'a')][0] / ref
axes[1].text(21.3, amp_aba * 1.2, r'$(a,b,a)$', color=GHT['a'], ha='left', va='center', fontsize=6.0)
axes[1].text(0.98, 0.94, 'uniform spacing', transform=axes[1].transAxes, ha='right', va='top', fontsize=6.6)

axes[1].set_xlabel('position along the fiber (m)')
axes[1].set_xlim(0, ZMAXPLOT)
fig.text(0.02, 0.55, 'correlation magnitude, normalized', rotation=90, ha='left', va='center', fontsize=7)

os.makedirs('figs', exist_ok=True)
fig.savefig('figs/fig_s51_ghostcorr.pdf')
fig.savefig('figs/fig_s51_ghostcorr.png', dpi=220)
print('saved figs/fig_s51_ghostcorr.pdf')
for let, Z in LAYOUTS.items():
    print('layout', let)
    for seq, amp, zpos in paths(Z):
        print('  %-10s amp %.2e  z %5.1f m' % ('(' + ','.join(seq) + ')', amp, zpos))
