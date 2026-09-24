"""Fig. 12 of the paper, in the layout of Fig. 14 of Markowski et al. (J. Lightwave Technol. 2023), from the decoded
coherent-model records of the 50-sensor array (jlt_chain.py, jlt_decode.py, cache/jlt/): top row the array as in [JLT] with a +-5-cm spacing tolerance, bottom row the designed array.
(a,c) reconstructed spectra of all 50 gratings, normalized to the peak of the first grating; (b,d) error of the restored Bragg wavelength
per grating in pm on a common scale (Gaussian fit over the +-S window around the nominal wavelength, reference-corrected in d) with the RMS error. The single-step correlograms of the earlier version are dropped (returns 0.75 chip apart are not resolved).
Output: figs/fig_jlt_style.pdf/.png"""
import sys, os
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
import figstyle as FS
import decode as D
FS.apply()
T = HERE / 'cache' / 'jlt'
C0 = 299792458.0
LAM_PM = 0.0        # sweep step shown in (a) and (d), offset from the nominal wavelength [pm]

def correlogram(rec, chip):
    """Correlogram at the sweep step nearest to the nominal wavelength (offset 0, the common Bragg wavelength of the
    uniform array and the center of the middle band), normalized to its largest peak, and its values at the
    nominal delays of the 50 gratings."""
    code = rec['code']; spc = int(rec['spc']); ng = float(rec['n_group']); x = rec['x_pm']
    m0 = int(np.argmin(np.abs(x - LAM_PM)))
    X = D.correlate(rec['det'][0, m0], code, spc)
    z_axis = np.arange(len(X)) / spc * C0 / (2 * ng * chip)          # sample -> equivalent position [m]
    bins = D.samples_of(rec['z'], ng, chip, spc); off = D.find_offset(np.array([X]), bins[0], spc)
    idx = ((np.round(bins).astype(int) + off) % len(X))
    top = X[idx].max()
    return z_axis, X / top, rec['z'], X[idx] / top


rows = [('J5_s3', 511 / 16.352e-6, 'uniform, co-tuned, $\\pm$5 cm', 'errc', 'Gaussian fit'),
        ('D4r_s1', 50e6, 'three bands, irregular', 'errc', 'Gaussian fit')]
fig, ax = plt.subplots(1, 3, figsize=(7.1, 1.85), layout='constrained', gridspec_kw=dict(width_ratios=[1.0, 1.0, 0.85]))
fig.get_layout_engine().set(wspace=0.10)
c = ax[2]
c.axhline(0, color='0.6', lw=0.5)
for r, (name, chip, lab, key, est) in enumerate(rows):
    dec = np.load(T / ('jltc_%s_dec.npz' % name))
    zs = dec['z']
    b = ax[r]
    S = np.clip(dec['S'], 0, None); x = dec['x_pm']
    im = b.pcolormesh(x / 1000.0, zs, S, cmap='viridis', vmin=0, vmax=1.3, shading='nearest', rasterized=True)
    b.set_xlim(x[0] / 1000.0, x[-1] / 1000.0); b.set_ylim(0, zs.max() + 2)
    b.set_xlabel('Wavelength offset [nm]')
    if r == 0:
        b.set_ylabel('Position [m]')
    else:
        b.set_yticklabels([])
    b.text(0.03, 0.95, lab, transform=b.transAxes, ha='left', va='top', fontsize=6, color='white')
    e = dec[key]
    k = np.arange(1, 51)
    col = FS.VERM if r == 0 else FS.BLUE
    c.plot(k, e, 'o' if r == 0 else 's', color=col, ms=2.4, mfc='none', mew=0.7, label='%s, RMS %.0f pm' % (lab.replace(', $\\pm$5 cm', ''), np.sqrt(np.mean(e ** 2))))
c.set_xlim(0, 51)
c.set_ylim(-150, 150)
c.set_ylabel('$\\lambda_B$ error [pm]')
c.set_xlabel('Grating number')
c.legend(fontsize=5.2, loc='upper left', handlelength=1.0, borderaxespad=0.3, labelspacing=0.2)
for a_, L in zip(ax, 'abc'):
    FS.letter(a_, L)
cb = fig.colorbar(im, ax=ax[:2], shrink=0.72, pad=0.035, aspect=16, anchor=(0.0, 0.0))
cb.ax.set_title('$\\dfrac{S_k}{\\max S_1}$', fontsize=6.5, pad=4, loc='center'); cb.ax.tick_params(labelsize=6)
OUT = HERE.parent / 'figs'
for ext in ('pdf', 'png'):
    fig.savefig(OUT / ('fig_jlt_style.' + ext), dpi=300)
print('saved fig_jlt_style')
