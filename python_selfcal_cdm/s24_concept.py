"""s24_concept.py - the whole system in one figure, v5.

Layout discipline modelled on Markowski et al., JLT 41(9) 2023, Fig. 1:
a compact block schematic with framed micro-oscillograms attached at their
tap points, and the signal-processing story in a tight row of panels below.

Two physical axes thread the figure: delay (= position) and wavelength.
  a  hardware chain: drive staircase (M steps, one code period each) ->
     swept VCSEL -> circulator -> K same-wavelength FBGs; echoes overlap on
     one photodiode; delay axis under the fiber states tau_k = 2 n z_k / c;
     frame strip gives the time cost and update rate.
  b  one sweep step: overlapping spectra, sampled by the step (dots).
  c  correlation of that step: peaks sorted by delay, heights equal to the
     sampled reflectivities of b, on the small code-leakage floor.
  d  M steps stack into the delay-wavelength map (full model with shadowing
     and leakage); row = step, column = grating.
  e  one column is one grating's sampled spectrum; a Gaussian fit returns
     lambda_B, whose shift is the measurand.

Colour continuity carries the argument: each grating keeps its colour from
the fiber, through the step, the correlation, the map, to the spectrum.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS

FS.apply()

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
NCH = 127
_MS = 1.0 - 2.0 * C._mls01(7)
_MS = np.roll(_MS, -int(np.argmax(_MS == 1.0)))
ACORR = C.periodic_xcorr(_MS, _MS)               # 1 at lag 0, -1/N elsewhere
ECOL = [FS.VERM, FS.ORANGE, FS.GREEN, FS.PURPLE]

TXT = '0.25'
GY = '0.45'

# --- model data --------------------------------------------------------------
M = 120
nu = np.linspace(-2.3 * F, 2.3 * F, M)
K = 4
bins = np.array([16, 47, 76, 105])
nub = np.array([-16.0, 9.0, -5.0, 19.0])
R = 0.05
shapes = np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
tcum = np.ones((K, M))
for k in range(1, K):
    tcum[k] = tcum[k - 1] * (1.0 - R * shapes[k - 1]) ** 2
prim = R * shapes * tcum
Wl = ACORR[(np.arange(NCH)[:, None] - bins[None, :]) % NCH]
grid = Wl @ prim                          # mainlobe + code leakage, full model
HL = 1                                    # highlighted grating (orange)
JM = np.argmin(np.abs(nu - 0.0))          # highlighted sweep step, lambda_m = 0
samp = shapes[:, JM]                      # unit-peak samples at that step
spec = grid[bins[HL]]
mu_fit = C.gauss_fit_peak(nu, spec)
lam = nu * PM / 1000.0
lam_m = nu[JM] * PM / 1000.0

# --- signals for the micro-insets -------------------------------------------
OSR = 6
taus = [3.0, 7.5, 12.0, 16.5]
drv2 = np.repeat(0.5 * (1 + _MS[:40]), OSR)
tt2 = np.arange(len(drv2)) / OSR
ech = []
for tcau, a in zip(taus, [1.0, 0.80, 0.62, 0.48]):
    e = np.zeros_like(drv2); sh = int(tcau * OSR)
    e[sh:] = a * drv2[:len(drv2) - sh]
    ech.append(e)
rngN = np.random.default_rng(3)
sm = np.sum(ech, axis=0) + rngN.normal(0, 0.06, len(drv2))
W2 = tt2 < 34

# staircase drive with one code period per step
NSTEP, CHPS = 4, 14
stair = []
for i in range(NSTEP):
    stair.append(i * 1.0 + 0.40 * 0.5 * (1 + _MS[i * CHPS:(i + 1) * CHPS]))
stair = np.repeat(np.concatenate(stair), OSR)
ts = np.arange(len(stair)) / (OSR * CHPS)

# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(7.1, 3.30))
bg = fig.add_axes([0, 0, 1, 1]); bg.set_xlim(0, 100); bg.set_ylim(0, 100)
FS.despine_all(bg)


def mini(rect, frame=True):
    a = fig.add_axes(rect)
    for sp in a.spines.values():
        sp.set_visible(frame)
        sp.set_linewidth(0.5)
        sp.set_color('0.55')
    a.set_xticks([]); a.set_yticks([])
    if frame:
        a.set_facecolor('white')
    return a


# =========================== a: the chain ===================================
YW = 87.0            # fiber height in bg units

# --- drive inset: staircase + one code period per step ---------------------
axdrv = mini([0.032, 0.735, 0.100, 0.170])
axdrv.plot(ts, stair, color=FS.BLUE, lw=0.6)
axdrv.set_xlim(-0.1, ts[-1] + 0.1); axdrv.set_ylim(-0.75, NSTEP + 0.15)
axdrv.annotate('', xy=(0.97, -0.35), xytext=(2.03, -0.35),
               arrowprops=dict(arrowstyle='<->', lw=0.5, color=GY))
axdrv.text(1.5, -0.62, 'code', fontsize=4.6, ha='center', va='top', color=GY)
axdrv.text(0.12, NSTEP - 0.05, '$I(t)$', fontsize=5.2, color=FS.BLUE,
           va='top')
bg.text(8.2, 97.6, 'drive: $M$ steps,', fontsize=5.4, ha='center', color=TXT)
bg.text(8.2, 94.4, 'one code period each', fontsize=5.4, ha='center',
        color=TXT)

# --- VCSEL box -------------------------------------------------------------
bg.annotate('', xy=(16.4, YW - 3), xytext=(13.6, YW - 3),
            arrowprops=dict(arrowstyle='-|>', lw=0.8, color=FS.BLUE))
bg.add_patch(FancyBboxPatch((16.8, YW - 7.2), 7.6, 8.4, fc='white',
                            ec='0.35', lw=0.7,
                            boxstyle='round,pad=0.3,rounding_size=0.9'))
bg.text(20.6, YW - 3, 'swept\nVCSEL', fontsize=5.6, ha='center', va='center',
        color='0.15')

# --- circulator ------------------------------------------------------------
CX = 28.4
bg.plot([24.8, CX - 1.7], [YW - 3, YW - 3], color='0.35', lw=0.9)
bg.plot([CX], [YW - 3], marker='o', ms=8.6, mfc='white', mec='0.35', mew=0.7)
th = np.linspace(0.6 * np.pi, 2.05 * np.pi, 40)
bg.plot(CX + 0.75 * np.cos(th), YW - 3 + 1.35 * np.sin(th), color='0.45',
        lw=0.5)
bg.plot([CX + 0.75 * np.cos(th[-1])], [YW - 3 + 1.35 * np.sin(th[-1])],
        marker='>', ms=1.4, color='0.45')

# --- fiber and gratings ----------------------------------------------------
YF = YW - 3
FX0, FX1 = CX + 1.7, 66.0
bg.plot([FX0, FX1], [YF, YF], color='0.35', lw=1.1, solid_capstyle='round')
bg.plot([FX1, FX1 + 3.2], [YF, YF], color='0.35', lw=1.1,
        ls=(0, (1.2, 1.6)))
gx = [34.5, 41.5, 47.5, 63.0]
for x, col, i in zip(gx, ECOL, range(4)):
    for dx in (-0.4, 0.0, 0.4):
        bg.plot([x + dx, x + dx], [YF - 2.3, YF + 2.3], color=col, lw=1.0)
    bg.text(x, YF + 3.2, '$z_%d$' % (i + 1), fontsize=5.4, ha='center',
            va='bottom', color=col)
bg.text(55.2, YF - 5.6, r'$\cdots$', fontsize=6.5, ha='center', va='bottom',
        color='0.35')
bg.text(71.0, YF + 3.0, '$K$ FBGs share one band:', fontsize=5.4, color=TXT,
        va='center')
bg.text(71.0, YF - 0.6, r'equal nominal $\lambda_B$,', fontsize=5.4,
        color=TXT, va='center')
bg.text(71.0, YF - 4.2, 'told apart by position', fontsize=5.4, color=TXT,
        va='center')

# --- delay axis under the fiber -------------------------------------------
YD = YF - 9.0
bg.annotate('', xy=(FX1 + 2.0, YD), xytext=(FX0, YD),
            arrowprops=dict(arrowstyle='-|>', lw=0.6, color=GY))
for x, col, i in zip(gx, ECOL, range(4)):
    bg.plot([x, x], [YF - 3.0, YD + 0.9], color=col, lw=0.45,
            ls=(0, (1.4, 1.4)))
    bg.plot([x, x], [YD - 0.8, YD + 0.8], color=col, lw=0.9)
    bg.text(x, YD - 2.0, r'$\tau_%d$' % (i + 1), fontsize=5.4, ha='center',
            va='top', color=col)
bg.text(71.0, YD - 1.6, r'round trip $\tau_k = 2 n z_k / c$:', fontsize=5.6,
        color=TXT, va='center')
bg.text(71.0, YD - 5.2, 'delay names the grating', fontsize=5.6, color=TXT,
        va='center')

# --- echo inset above the fiber -------------------------------------------
axech = mini([0.500, 0.870, 0.120, 0.115])
for e, col, off in zip(ech, ECOL, (3.6, 2.4, 1.2, 0.0)):
    axech.plot(tt2[W2], off + e[W2], color=col, lw=0.55)
axech.set_xlim(0, 34); axech.set_ylim(-0.4, 5.1)
bg.text(56.0, 99.5, 'echoes: same code, delayed, weighted', fontsize=5.4,
        color=TXT, ha='center')

# --- photodiode, sum, correlator -------------------------------------------
bg.plot([CX, CX], [YF - 1.9, YF - 10.6], color='0.35', lw=0.8)
bg.plot([CX], [YF - 12.1], marker='v', ms=4.0, color='0.35')
bg.plot([CX - 1.0, CX + 1.0], [YF - 14.2, YF - 14.2], color='0.35', lw=0.8)
bg.text(CX + 1.9, YF - 12.4, 'PD', fontsize=5.6, color=TXT, ha='left',
        va='center')
axsum = mini([0.085, 0.525, 0.100, 0.115])
axsum.plot(tt2[W2], sm[W2], color='0.3', lw=0.55)
axsum.set_xlim(0, 34)
bg.text(13.5, 65.5, 'one photodiode sums them', fontsize=5.4, color=TXT,
        ha='center')
bg.plot([18.7, 26.9], [61.0, YF - 12.6], color='0.65', lw=0.45)

bg.plot([CX, CX], [YF - 15.9, YF - 18.6], color='0.35', lw=0.8)
bg.add_patch(FancyBboxPatch((CX - 5.5, YF - 24.0), 11.0, 5.4, fc='white',
                            ec='0.35', lw=0.7,
                            boxstyle='round,pad=0.25,rounding_size=0.8'))
bg.text(CX, YF - 21.3, 'correlate with $c$', fontsize=5.6, ha='center',
        va='center', color='0.15')

# --- frame strip: the time cost of one measurement -------------------------
FSX0, FSX1, FSY0, FSY1 = 45.5, 70.0, 59.0, 63.6
bg.add_patch(plt.Rectangle((FSX0, FSY0), FSX1 - FSX0, FSY1 - FSY0,
                           fc='white', ec='0.45', lw=0.6))
cw = (FSX1 - FSX0) / 5.0
for i in range(1, 5):
    bg.plot([FSX0 + i * cw, FSX0 + i * cw], [FSY0, FSY1], color='0.45',
            lw=0.5)
bg.add_patch(plt.Rectangle((FSX0, FSY0), cw, FSY1 - FSY0, fc=FS.BLUE,
                           alpha=0.18, lw=0))
for i, lab in enumerate(('1', '2', '3', r'$\cdots$', '$M$')):
    bg.text(FSX0 + (i + 0.5) * cw, 0.5 * (FSY0 + FSY1), lab, fontsize=5.2,
            ha='center', va='center', color='0.3')
bg.text(0.5 * (FSX0 + FSX1), FSY1 + 1.8,
        r'one frame: $M$ steps $\times$ $N$ chips', fontsize=5.4,
        ha='center', color=TXT)
bg.text(0.5 * (FSX0 + FSX1), FSY0 - 2.0,
        r'update rate $f_r = f_\mathrm{chip}/(MN)$', fontsize=5.6,
        ha='center', va='top', color=TXT)

# --- colour key ------------------------------------------------------------
for i, col in enumerate(ECOL):
    bg.plot([76.5 + 1.9 * i], [62.6], marker='o', ms=2.8, color=col)
bg.text(84.5, 62.6, r'$=$ grating identity,', fontsize=5.2, ha='left',
        va='center', color=TXT)
bg.text(76.5, 59.0, 'kept through every panel', fontsize=5.2, ha='left',
        va='center', color=TXT)

FS.panel(fig, 0.006, 0.995, 'a')

# =========================== bottom row =====================================
BY, BH = 0.112, 0.305
axb = fig.add_axes([0.055, BY, 0.160, BH])
axc = fig.add_axes([0.293, BY, 0.180, BH])
axd = fig.add_axes([0.553, BY, 0.210, BH])
axe = fig.add_axes([0.855, BY, 0.132, BH])

# --- b: one sweep step samples every grating -------------------------------
for c0, col, s in zip(nub, ECOL, samp):
    g = np.exp(-0.5 * ((lam * 1000 / PM - c0) / SIG) ** 2)
    axb.plot(lam, g, color=col, lw=0.8)
    axb.plot([lam_m], [s], marker='o', ms=2.6, color=col, zorder=5)
axb.axvline(lam_m, color='0.2', lw=0.6, ls=(0, (3, 1.5)))
axb.text(lam_m + 0.035, 1.16, r'step $\lambda_m$', fontsize=5.4,
         color='0.15')
axb.set_xlim(-0.50, 0.50); axb.set_ylim(0, 1.30)
axb.set_xticks([-0.4, 0, 0.4]); axb.set_yticks([0, 1])
axb.tick_params(labelsize=5.4, pad=1.5)
axb.set_xlabel(r'$\lambda$ offset (nm)', fontsize=5.8, labelpad=1.5)
axb.set_ylabel(r'reflectivity $R_k(\lambda)$', fontsize=5.8, labelpad=1)

# --- c: the correlation of that step (unshadowed, to match the dots of b) ---
tr = (Wl @ (R * shapes[:, JM])) / R
axc.fill_between(np.arange(NCH), 0, np.abs(tr), color='0.78', lw=0)
for b, col in zip(bins, ECOL):
    h = tr[b]
    axc.plot([b, b], [0, h], color=col, lw=1.4, solid_capstyle='butt')
    axc.plot([b], [h], marker='o', ms=2.6, color=col)
axc.annotate('height $= R_k(\\lambda_m)$',
             xy=(bins[2] + 1.5, tr[bins[2]]), xytext=(60, 1.12),
             fontsize=5.4, color='0.15',
             arrowprops=dict(arrowstyle='-', lw=0.5, color='0.6'))
axc.annotate('code leakage,\n' + r'$\sim\!1/N$',
             xy=(92, 0.035), xytext=(97, 0.44), fontsize=5.0, color=GY,
             ha='center',
             arrowprops=dict(arrowstyle='-', lw=0.5, color='0.6'))
axc.set_xlim(0, NCH); axc.set_ylim(0, 1.30)
axc.set_xticks(list(bins))
axc.set_xticklabels([r'$\tau_1$', r'$\tau_2$', r'$\tau_3$', r'$\tau_4$'],
                    fontsize=5.6)
for t, col in zip(axc.get_xticklabels(), ECOL):
    t.set_color(col)
axc.set_yticks([0, 1])
axc.tick_params(labelsize=5.4, pad=1.5)
axc.set_xlabel(r'delay $\tau$', fontsize=5.8, labelpad=1.5)
axc.set_ylabel('correlation', fontsize=5.8, labelpad=1)

# --- d: the map -------------------------------------------------------------
axd.imshow(grid.T, aspect='auto', origin='lower', cmap='Greys',
           extent=[0, NCH, lam[0], lam[-1]], vmax=grid.max() * 0.9)
for b, col in zip(bins, ECOL):
    axd.plot(b, lam[-1] + 0.022, marker='v', ms=2.8, color=col,
             clip_on=False)
axd.axvspan(bins[HL] - 3.5, bins[HL] + 3.5, color=FS.ORANGE, alpha=0.16,
            lw=0)
axd.axhline(lam_m, color='0.2', lw=0.55, ls=(0, (3, 1.5)))
axd.text(3.0, lam_m + 0.035, r'row $=$ step $\lambda_m$', fontsize=5.0,
         color='0.15')
axd.text(bins[HL] + 6.0, lam[0] + 0.03, 'column $=$ grating', fontsize=5.0,
         color='0.15', ha='left', va='bottom')
axd.set_xlim(0, NCH)
axd.set_xticks(list(bins))
axd.set_xticklabels([r'$\tau_1$', r'$\tau_2$', r'$\tau_3$', r'$\tau_4$'],
                    fontsize=5.6)
for t, col in zip(axd.get_xticklabels(), ECOL):
    t.set_color(col)
axd.set_yticks([-0.4, 0, 0.4])
axd.tick_params(labelsize=5.4, pad=1.5)
axd.set_xlabel(r'delay $\tau$', fontsize=5.8, labelpad=1.5)
axd.set_ylabel(r'$\lambda$ offset (nm)', fontsize=5.8, labelpad=0.5)
for sp in ('top', 'right'):
    axd.spines[sp].set_visible(False)

# --- e: the readout ---------------------------------------------------------
axe.plot(lam, spec / spec.max(), color=ECOL[HL], lw=0.9)
axe.plot(lam[::6], (spec / spec.max())[::6], ls='none', marker='o', ms=1.7,
         color=ECOL[HL])
axe.axvline(mu_fit * PM / 1000, color='0.35', ls=(0, (3, 2)), lw=0.7)
axe.text(mu_fit * PM / 1000 + 0.06, 0.80, r'$\lambda_B$', fontsize=6.4,
         color='0.15', va='top')
axe.text(0.05, 0.97, r'$\Delta\lambda_B \propto \Delta T,\ \varepsilon$',
         fontsize=5.8, color='0.15', transform=axe.transAxes, va='top')
axe.set_xlim(-0.50, 0.50); axe.set_ylim(0, 1.30)
axe.set_xticks([-0.4, 0, 0.4]); axe.set_yticks([0, 1])
axe.tick_params(labelsize=5.4, pad=1.5)
axe.set_xlabel(r'$\lambda$ offset (nm)', fontsize=5.8, labelpad=1.5)
axe.set_ylabel('refl.', fontsize=5.8, labelpad=1)

# --- narrative titles over the bottom row ----------------------------------
TY = BY + BH + 0.048
fig.text(0.135, TY, 'one step samples every grating', fontsize=5.6,
         ha='center', va='center', color='0.15')
fig.text(0.383, TY, 'correlation sorts by delay', fontsize=5.6,
         ha='center', va='center', color='0.15')
fig.text(0.658, TY, r'$M$ steps $\rightarrow$ delay$-\lambda$ map',
         fontsize=5.6, ha='center', va='center', color='0.15')
fig.text(0.921, TY, r'peak $\rightarrow$ measurand', fontsize=5.6,
         ha='center', va='center', color='0.15')

# --- flow arrows ------------------------------------------------------------
bg.annotate('', xy=(34.5, 51.5), xytext=(CX + 6.0, YF - 21.3),
            arrowprops=dict(arrowstyle='-|>', lw=1.0, color=GY,
                            connectionstyle='arc3,rad=-0.25'))
bg.text(38.3, 56.5, 'one step', fontsize=5.2, color=GY, ha='center')
YA = (BY + 0.5 * BH) * 100
bg.annotate('', xy=(26.6, YA), xytext=(22.2, YA),
            arrowprops=dict(arrowstyle='-|>', lw=0.9, color=GY))
bg.annotate('', xy=(52.0, YA), xytext=(48.0, YA),
            arrowprops=dict(arrowstyle='-|>', lw=0.9, color=GY))
for i in range(3):
    bg.plot([48.3 + 0.4 * i, 51.0 + 0.4 * i],
            [YA + 3.4 + 1.1 * i, YA + 3.4 + 1.1 * i], color='0.7', lw=0.55)
bg.text(50.0, YA + 8.2, r'$\times M$', fontsize=5.2, color=GY, ha='center')
bg.annotate('', xy=(83.6, YA), xytext=(79.4, YA),
            arrowprops=dict(arrowstyle='-|>', lw=0.9, color=GY))

PLY = TY + 0.045
FS.panel(fig, 0.010, PLY, 'b')
FS.panel(fig, 0.248, PLY, 'c')
FS.panel(fig, 0.508, PLY, 'd')
FS.panel(fig, 0.812, PLY, 'e')

fig.savefig('figs/fig_s24_concept.png', bbox_inches='tight', pad_inches=0.02)
print('fitted lambda_B of the highlighted grating: %.1f pm (true %.1f pm)'
      % (mu_fit * PM, nub[HL] * PM))
print('saved figs/fig_s24_concept.png')
