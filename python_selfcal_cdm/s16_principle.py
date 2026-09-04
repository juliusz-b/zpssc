"""s16_principle.py - the two explanatory figures of the paper.

Figure 1, principle: the optical layout on the left and, on the right, what the
correlator actually returns. The despread output is a two-dimensional map, delay
bin against wavelength: the code marks WHICH wavelength was launched, the delay
marks WHICH grating sent the light back. A horizontal cut through that map is one
grating's reflection spectrum, and its peak is the measurement.

Figure 2, mechanisms: the three array effects that this study quantifies, each
drawn from the same model used for the results, not sketched by hand.
  (a) multiple reflections: a third-order path returns at tau_a - tau_b + tau_c,
      which for uniform spacing is always another occupied bin and for randomised
      spacing usually is not;
  (b) code leakage: the autocorrelation side lobe adds a scaled copy of every
      other grating's spectrum as a broad background under the wanted line.

Spectral shadowing, the third array term, has its own figure in s19 because it
comes with a correction.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS
FS.apply()

PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
NCH = 127
MSEQ = 1.0 - 2.0 * C._mls01(7)
ACORR = C.periodic_xcorr(MSEQ, MSEQ)


# ===========================================================================
# FIGURE 1: principle
# ===========================================================================
fig = plt.figure(figsize=(7.1, 2.65))
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.28)

# --- (a) optical layout ----------------------------------------------------
ax = fig.add_subplot(gs[0, 0]); ax.axis('off')
ax.set_xlim(0, 10.4); ax.set_ylim(0, 6.2)


def box(x, y, w, h, t, fc='#eaf1fb', fs=7.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle='round,pad=0.05,rounding_size=0.12',
                                fc=fc, ec='#0072B2', lw=1.2))
    ax.text(x + w / 2, y + h / 2, t, ha='center', va='center', fontsize=fs)


def arrow(x1, y1, x2, y2, col='#333', lw=1.2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>',
                                 mutation_scale=10, lw=lw, color=col))


box(0.15, 4.20, 2.55, 1.35, 'swept VCSEL,\ncode-modulated', fs=7.4)
# code waveform inside the source box area
tt = np.linspace(0.35, 2.35, 300)
bits = MSEQ[:16]
cw = np.repeat(bits, len(tt) // 16 + 1)[:len(tt)]
ax.plot(tt, 3.55 + 0.18 * cw, color='#0072B2', lw=0.9)
ax.text(1.35, 3.15, 'code c(t)', ha='center', fontsize=7, color='#0072B2')

ax.add_patch(Circle((3.5, 4.85), 0.42, fc='#fdf3e7', ec='#0072B2', lw=1.2))
ax.text(3.5, 4.85, 'circ', ha='center', va='center', fontsize=7)
arrow(2.72, 4.85, 3.05, 4.85)

# fiber with gratings at irregular positions
ax.plot([3.95, 10.2], [4.85, 4.85], color='#555', lw=1.6)
zpos = [4.7, 5.9, 7.6, 8.5, 9.7]
for i, z in enumerate(zpos):
    ax.add_patch(Rectangle((z - 0.16, 4.62), 0.32, 0.46, fc='#D55E00',
                           ec='#7b241c', hatch='///', lw=0.8))
    ax.text(z, 5.30, r'$z_%d$' % (i + 1), ha='center', fontsize=7)
ax.annotate('', xy=(4.7, 4.35), xytext=(5.9, 4.35),
            arrowprops=dict(arrowstyle='<->', lw=0.8, color='#555'))
ax.annotate('', xy=(7.6, 4.35), xytext=(8.5, 4.35),
            arrowprops=dict(arrowstyle='<->', lw=0.8, color='#555'))
ax.text(7.2, 3.95, 'randomized spacing (Section IV-B)', ha='center', fontsize=7,
        color='#555')
ax.text(7.1, 5.75, 'FBG array: same nominal wavelength,\nseparated only by delay',
        ha='center', fontsize=7.4)

# return path
arrow(3.5, 4.43, 3.5, 3.0)
box(2.75, 2.05, 1.5, 0.9, 'InGaAs\nAPD')
arrow(4.25, 2.5, 5.05, 2.5)
box(5.00, 2.05, 2.45, 0.9, 'equivalent-time\nsampling', fs=7.2)
arrow(7.45, 2.5, 8.15, 2.5)
box(8.15, 2.05, 2.2, 0.9, 'matched-filter\ncorrelator', fs=7.2)
arrow(9.25, 2.05, 9.25, 1.35)
box(6.6, 0.35, 3.7, 1.0, 'delay-wavelength map  ->  peak fit  ->  $\\lambda_B$',
    fc='#eafaf0')

ax.text(0.1, 1.3, 'one code on the fiber\nat a time: the sweep\nvisits wavelengths\n'
                  'in sequence, so gratings\nare separated by the\nAUTOCORRELATION\n'
                  'side lobe', fontsize=7.2, va='center', color='#0072B2')
FS.letter(ax, 'a')
# --- (b) despread map ------------------------------------------------------
ax2 = fig.add_subplot(gs[0, 1])
M = 96
nu = np.linspace(-2.6 * F, 2.6 * F, M)
rng = np.random.default_rng(2)
K = 5
bins = np.array([11, 27, 46, 79, 103])
nub = np.array([-18.0, 9.0, -4.0, 21.0, -11.0])
R = 0.05
shapes = np.exp(-0.5 * ((nu[None, :] - nub[:, None]) / SIG) ** 2)
tcum = np.ones((K, M))
for k in range(1, K):
    tcum[k] = tcum[k - 1] * (1.0 - R * shapes[k - 1]) ** 2
prim = R * shapes * tcum
grid = np.zeros((NCH, M))
grid[bins] = prim
W = ACORR[(np.arange(NCH)[:, None] - bins[None, :]) % NCH]
grid = grid + W @ prim
im = ax2.imshow(grid.T, aspect='auto', origin='lower', cmap='viridis',
                extent=[0, NCH, nu[0] * PM / 1000.0, nu[-1] * PM / 1000.0])
for k in range(K):
    ax2.plot(bins[k], nub[k] * PM / 1000.0, 'o', mfc='none', mec='w', ms=9, mew=1.2)
    ax2.text(bins[k] + 2, nub[k] * PM / 1000.0 + 0.06, 'FBG %d' % (k + 1),
             color='w', fontsize=7)
ax2.set_xlabel('delay bin  [chip]  ->  grating position')
ax2.set_ylabel('wavelength offset [nm]')
FS.letter(ax2, 'b')
cb = fig.colorbar(im, ax=ax2, pad=0.02); cb.set_label('despread reflectance', fontsize=7.5)
cb.ax.tick_params(labelsize=7)
ax2.text(0.98, 0.04, 'a horizontal cut is one grating spectrum',
         transform=ax2.transAxes, ha='right', fontsize=7, color='w')
fig.savefig('figs/fig_s16_principle.png', dpi=150, bbox_inches='tight')

# ===========================================================================
# FIGURE 2: error mechanisms
# ===========================================================================
# This figure is included at 0.72 text width in the paper. Draw it at that
# physical width so the 7 pt labels remain 7 pt after LaTeX placement.
fig2, ax = plt.subplots(1, 2, figsize=(5.1, 1.98))

# --- (b) ghost delays ------------------------------------------------------
axb = ax[0]
axb.set_xlim(-2, 102); axb.set_ylim(-0.35, 10.9); axb.axis('off')
DOT = 0.24                          # odstep kropek w slupku

UNI = [10, 25, 40, 55, 70]          # skok 15
RND = [15, 25, 41, 67, 86]          # odstepy 10, 16, 26, 19


def in_span_ghosts(bins):
    """Wszystkie sciezki trzeciego rzedu, ktore wracaja w obrebie tablicy."""
    out = []
    for ib, b in enumerate(bins):
        for a in bins[ib + 1:]:
            for c in bins[ib + 1:]:
                g = a - b + c
                if bins[0] <= g <= bins[-1]:
                    out.append(g)
    return out


def draw_row(y, bins, col, title):
    paths = in_span_ghosts(bins)
    mult = Counter(paths)
    hits = [x for x in paths if x in bins]
    # wlokno z siatkami: fizyczny uklad, ten sam w obu wierszach
    axb.plot([0, 100], [y, y], color='#555', lw=1.2, zorder=1)
    for x in bins:
        axb.add_patch(Rectangle((x - 1.5, y - 0.26), 3.0, 0.52, fc=col,
                                ec='#333', lw=0.5, zorder=3))
    # jedna kropka na sciezke, slupek na wspolne opoznienie
    for x in sorted(mult):
        on = x in bins
        c = '#D55E00' if on else '0.62'
        axb.plot([x, x], [y - 0.30, y - 0.70], color=c,
                 lw=1.0 if on else 0.7, zorder=2)
        for i in range(mult[x]):
            axb.plot([x], [y - 0.84 - DOT * i], marker='o', ms=2.3,
                     color=c, mec='none', zorder=4)
    axb.text(0, y + 0.62, title, fontsize=6.8, color=col)
    # podpis pod najglebszym slupkiem wiersza, bo slupki maja rozna dlugosc
    deep = y - 0.84 - DOT * (max(mult.values()) - 1)
    axb.text(0, deep - 0.78, '%d ghost paths return inside the array, %d on a grating'
             % (len(paths), len(hits)), fontsize=6.2,
             color='#D55E00' if hits else '0.45')


axb.text(0, 10.45, r'a ghost returns at $\tau_a-\tau_b+\tau_c$',
         fontsize=6.4, color='0.35')
draw_row(9.00, UNI, '#D55E00', 'uniform spacing: every ghost lands on a grating')
draw_row(4.35, RND, '#0072B2', 'randomized spacing: almost none does')

axb.annotate('', xy=(100, 0.85), xytext=(0, 0.85),
             arrowprops=dict(arrowstyle='-|>', color='0.55', lw=0.7))
axb.text(50, 0.02, 'delay bin $\\tau$ (grating position)', ha='center',
         fontsize=6.8, color='0.4')
FS.letter(axb, 'a')
# --- (c) code leakage ------------------------------------------------------
Mc = 96
nuc = np.linspace(-2.6 * F, 2.6 * F, Mc)
rng3 = np.random.default_rng(7)
Kc = 32
nubs = rng3.uniform(-25, 25, Kc)
Ac = 0.05 * np.exp(-0.5 * ((nuc[None, :] - nubs[:, None]) / SIG) ** 2)
wanted = Ac[0]
# real leakage for both codes: every grating in its own random delay bin,
# and the periodic autocorrelation of the transmitted sequence at the bin
# difference weights the spectrum of every other grating (eq. leakterm)
binsc = np.sort(rng3.choice(np.arange(1, NCH), size=Kc, replace=False))


def autocorr(code01):
    code01 = np.asarray(code01)
    pm = code01.astype(float) if code01.min() < 0 else C._to_pm1(code01).astype(float)   # gold_set is already bipolar
    r = np.fft.ifft(np.fft.fft(pm) * np.conj(np.fft.fft(pm))).real
    return r / r[0]


rho_m = autocorr(C._mls01(7))                 # -1/N at every nonzero lag
rho_g = autocorr(C.gold_set(7)[2])            # three values of mixed sign
print('Gold side lobes at N=%d: %s' % (NCH, sorted(set(np.round(rho_g[1:] * NCH).astype(int)))))


def leakage(rho):
    out = np.zeros_like(wanted)
    for j in range(1, Kc):
        out += rho[(binsc[0] - binsc[j]) % NCH] * Ac[j]
    return out


leak_m, leak_g = leakage(rho_m), leakage(rho_g)
xnm = nuc * PM / 1000.0
c_true = C.gauss_fit_peak(nuc, wanted) * PM
c_m = C.gauss_fit_peak(nuc, wanted + leak_m) * PM
c_g = C.gauss_fit_peak(nuc, wanted + leak_g) * PM
print('leakage panel: fitted center moves by %.1f pm (m-sequence) and %.1f pm (Gold), K=%d, N=%d'
      % (c_m - c_true, c_g - c_true, Kc, NCH))
ax[1].axhline(0, color='0.75', lw=0.6, zorder=1)
ax[1].plot(xnm, wanted / 0.05, color='#0072B2', lw=1.5, label='$A_k$, the line of grating $k$')
ax[1].plot(xnm, leak_m / 0.05, color='0.45', lw=1.1, label='$L_k$, m-sequence')
ax[1].plot(xnm, (wanted + leak_m) / 0.05, color='0.2', lw=1.0, ls=(0, (3, 1.5)), label='$A_k+L_k$, m-sequence')
ax[1].plot(xnm, leak_g / 0.05, color='#D55E00', lw=1.1, label='$L_k$, Gold code')
ax[1].plot(xnm, (wanted + leak_g) / 0.05, color='#D55E00', lw=1.0, ls=(0, (3, 1.5)), label='$A_k+L_k$, Gold code')
for cc, col in ((c_true, '#0072B2'), (c_m, '0.2'), (c_g, '#D55E00')):
    ax[1].plot([cc / 1000.0, cc / 1000.0], [1.04, 1.16], color=col, lw=0.9)
ax[1].set_ylim(-0.6, 2.15)
ax[1].set_yticks([-0.5, 0, 0.5, 1.0])
ax[1].set_xlabel('wavelength offset [nm]'); ax[1].set_ylabel('reflectance / $R$')
FS.letter(ax[1], 'b')
ax[1].grid(False, which='both', alpha=0.2)
ax[1].legend(fontsize=5.2, loc='upper left', ncol=2, frameon=True, columnspacing=0.8,
             handlelength=1.6, labelspacing=0.2, borderaxespad=0.3)

fig2.subplots_adjust(left=0.055, right=0.99, top=0.87, bottom=0.19,
                     wspace=0.30)
fig2.savefig('figs/fig_s16_mechanisms.png', dpi=150, bbox_inches='tight')
fig2.savefig('figs/fig_s16_mechanisms.pdf', bbox_inches='tight')

# --- kontrola liczbowa: statystyka zbiorcza cytowana w III-B --------------
# Rysunek pokazuje jednego ducha, bo to on niesie mechanizm. Skala zjawiska
# jest w tekscie, wiec liczona jest tutaj, na dziesieciu siatkach.
Kg = 10
bu = 3 + 7 * np.arange(Kg)
br = np.sort(np.random.default_rng(11).choice(np.arange(1, 90), size=Kg,
                                              replace=False))


def ghost_bins(b):
    out = []
    for a in range(len(b)):
        for bb in range(len(b)):
            for c in range(len(b)):
                if bb < a and bb < c:
                    gb = b[a] - b[bb] + b[c]
                    if b.min() <= gb <= b.max():
                        out.append(gb)
    return np.array(out)


print('ghost bins, uniform: %d in span, %.0f%% on a grating bin'
      % (len(ghost_bins(bu)), 100 * np.mean([g in bu for g in ghost_bins(bu)])))
print('ghost bins, random : %d in span, %.0f%% on a grating bin'
      % (len(ghost_bins(br)), 100 * np.mean([g in br for g in ghost_bins(br)])))
print('saved figs/fig_s16_principle.png, figs/fig_s16_mechanisms.png')
