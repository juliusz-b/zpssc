"""s21_cdmwdm.py - how the hybrid CDM-WDM system works, and why banding helps.

The largest coded arrays in the literature are hybrids: wavelength division
splits the array into bands, and code/delay separates the gratings inside each
band. This figure explains the addressing and then checks, numerically, a point
that the demonstration papers do not spell out: banding does not merely add
capacity, it RESETS the error mechanisms at every band edge.

  (a) The addressing plane. Every grating owns a coordinate pair: its wavelength
      band and its delay bin. Two gratings may share a delay bin as long as they
      sit in different bands, and share a band as long as they sit in different
      bins. The delay axis is therefore reused once per band.

  (b) The measurement cycle in time. The sweep is a staircase over the
      wavelength samples; at each step the code runs for the few sequence
      periods that equivalent-time sampling needs; the correlator then owns one
      column of the map. The frame time is M steps times the per-step
      acquisition, which sets the refresh rate of every sensor at once.

  (c) The numerical check. All three array mechanisms need spectral overlap:
      shadowing is a notch cut at the upstream grating's wavelength, a
      third-order ghost needs all three gratings to overlap, and leakage only
      pulls the fit when the leaked line lands near the fitted one. Gratings in
      another band are spectrally distant, so none of the three couples across
      bands. Consequence: an array of 32 gratings split into 4 bands of 8 should
      behave like a single band of 8, not like a single band of 32. The panel
      verifies exactly that, with the full model (shadowing, ghosts, leakage,
      noise).
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import warnings; warnings.filterwarnings('ignore')
import common as C
import figstyle as FS
FS.apply()


PM = C.PM_PER_GHZ
F = C.FBG_FWHM_GHZ
SIG = F / 2.35482
NCH = 127
_MS = 1.0 - 2.0 * C._mls01(7)
ACORR = C.periodic_xcorr(_MS, _MS)
SIGMA_N = 1.1e-6
DETUNE = 25.0                     # +-200 pm detuning inside a band
BAND_STEP = 8.0 * F               # band pitch, 2 nm: no spectral interaction

BAND_COLS = ['#0072B2', '#009E73', '#E69F00', '#CC79A7']


# ---------------------------------------------------------------------------
# full array model with bands (shadowing, ghosts, leakage, noise)
# ---------------------------------------------------------------------------
def run_banded(K, W, R, rng):
    """RMS error over K gratings split into W wavelength bands.

    A wide grid covers all bands; shadowing, ghosts and leakage all act on it,
    so any cross-band coupling the physics allows is present in the model."""
    per = K // W
    centres = BAND_STEP * (np.arange(W) - (W - 1) / 2.0)
    nub = np.concatenate([c + rng.uniform(-DETUNE, DETUNE, per) for c in centres])
    order = rng.permutation(K)               # physical order along the fiber
    nub = nub[order]
    band = np.repeat(np.arange(W), per)[order]
    span = BAND_STEP * W / 2.0 + 3 * F
    Mw = int(24 * span / F)                  # keep ~4 samples per sigma
    nuw = np.linspace(-span, span, Mw)
    shapes = np.exp(-0.5 * ((nuw[None, :] - nub[:, None]) / SIG) ** 2)
    refl = R * shapes
    tcum = np.ones((K, Mw))
    for k in range(1, K):
        tcum[k] = tcum[k - 1] * (1.0 - refl[k - 1]) ** 2
    primary = refl * tcum
    bins = np.sort(rng.choice(np.arange(1, NCH), size=K, replace=False))
    # third-order ghosts on occupied bins
    idx = np.arange(K)
    a, b, c = np.meshgrid(idx, idx, idx, indexing='ij')
    m = (b < a) & (b < c)
    a, b, c = a[m], b[m], c[m]
    gbin = bins[a] - bins[b] + bins[c]
    pos = np.clip(np.searchsorted(bins, gbin), 0, K - 1)
    hit = bins[pos] == gbin
    a, b, c, tgt = a[hit], b[hit], c[hit], pos[hit]
    cbar = (nub[a] + nub[b] + nub[c]) / 3.0
    spread = (nub[a] - cbar) ** 2 + (nub[b] - cbar) ** 2 + (nub[c] - cbar) ** 2
    amp = R ** 3 * np.exp(-0.5 * spread / SIG ** 2)
    G = np.zeros((K, Mw))
    sig_g = SIG / np.sqrt(3.0)
    for i in range(len(amp)):
        if amp[i] > 1e-12:
            G[tgt[i]] += amp[i] * np.exp(-0.5 * ((nuw - cbar[i]) / sig_g) ** 2)
    Wl = ACORR[(bins[:, None] - bins[None, :]) % NCH]
    np.fill_diagonal(Wl, 0.0)
    S = primary + Wl @ primary + G + rng.normal(0, SIGMA_N, (K, Mw))
    errs = np.empty(K)
    for k in range(K):
        w = np.abs(nuw - nub[k]) < 3.5 * F          # the sweep window of its band
        errs[k] = (C.gauss_fit_peak(nuw[w], S[k][w]) - nub[k]) * PM
    return float(np.sqrt(np.mean(errs ** 2)))


def avg(K, W, R, seed, ntrials=6):
    return float(np.mean([run_banded(K, W, R, np.random.default_rng(seed + t))
                          for t in range(ntrials)]))


R0 = 0.10
e_one32 = avg(32, 1, R0, 100)
e_band = avg(32, 4, R0, 200)
e_one8 = avg(8, 1, R0, 300)

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(4.8, 2.20))

# --- (a) addressing plane ---------------------------------------------------
axa = ax[0]
rng = np.random.default_rng(4)
Wb = 4
for wί in range(Wb):
    axa.axhspan(wί, wί + 1, color=BAND_COLS[wί], alpha=0.10)
    axa.text(2, wί + 0.82, 'band %d' % (wί + 1), fontsize=6.8,
             color=BAND_COLS[wί])
SHARED = 84                      # one delay bin deliberately reused by every band
pool = np.array([x for x in range(5, NCH - 5) if abs(x - SHARED) > 6])
for wί in range(Wb):
    b = np.sort(np.append(rng.choice(pool, size=6, replace=False), SHARED))
    lam = rng.uniform(0.15, 0.72, size=b.size)
    axa.plot(b, wί + lam, 'v', ms=5, color=BAND_COLS[wί])
axa.axvline(SHARED, color='0.3', ls=':', lw=0.9)
axa.annotate('one bin, four bands', xy=(SHARED, 3.92), xytext=(SHARED - 6, 4.16),
             fontsize=5.8, color='0.3', ha='right',
             arrowprops=dict(arrowstyle='-', color='0.3', lw=0.5))
axa.set_xlim(0, NCH); axa.set_ylim(0, Wb + 0.30)
axa.set_yticks([])
axa.set_xlabel('delay bin (grating position)')
axa.set_ylabel('wavelength band')
axa.set_title('(a) CDM-WDM addressing', fontsize=9)

# --- (b) banding resets the mechanisms ---------------------------------------
axc = ax[1]
vals = [e_one32, e_band, e_one8]
labels = ['one band' '\n' r'$K = 32$',
          '4 bands' '\n' 'of 8' '\n' r'($K = 32$)',
          'one band' '\n' r'$K = 8$' '\n' '(reference)']
cols = ['#D55E00', '#E69F00', '#0072B2']
bars = axc.bar(np.arange(3), vals, 0.6, color=cols)
for i, v in enumerate(vals):
    axc.text(i, v * 1.04, '%.1f pm' % v, ha='center', fontsize=7.5)
axc.set_xticks(np.arange(3)); axc.set_xticklabels(labels, fontsize=7)
axc.set_ylabel('RMS Bragg error [pm]')
axc.set_ylim(0, max(vals) * 1.22)
axc.set_title('(b) Banding resets the error chain', fontsize=9)
axc.set_axisbelow(True)      # siatka pod slupkami, nie na nich
axc.grid(True, axis='y', alpha=0.25)

fig.tight_layout()
fig.savefig('figs/fig_s21_cdmwdm.png', dpi=150, bbox_inches='tight')
fig.savefig('figs/fig_s21_cdmwdm.pdf', bbox_inches='tight')

print('full model, R = %.0f%%: one band K=32: %.2f pm; 4 bands x 8: %.2f pm;'
      ' one band K=8: %.2f pm' % (R0 * 100, e_one32, e_band, e_one8))
print('=> the banded array behaves like its band, not like its total:')
print('   shadowing, ghosts and leakage all need spectral overlap, so none of')
print('   them couples across bands; capacity multiplies by the band count')
print('saved figs/fig_s21_cdmwdm.png')
