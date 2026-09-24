"""50-sensor array of [JLT] (Markowski et al., J. Lightwave Technol. 41, 2892, 2023) in the coherent model,
variants with source realizations and references. Records go to records/jlt/, decode them with jlt_decode.py.
Usage: python jlt_chain.py <variant> <seed>
Variants: D0 (as in [JLT]: co-tuned, uniform 2.5 m, 31.25 Mchip/s, random code), D1 (+ m-sequence), D2 (+ 50 Mchip/s),
D3 (+ irregular positions), D4 (+ three 0.8-nm bands), D4r (D4 + two references at 128 and 131 m on a separate branch,
Bragg wavelengths at -950 and +950 pm), D3r, D0r (references added to D3 / D0 for the same correction).
The quick decoding printed here (local Gaussian fit) is only a progress check. The estimator of the paper is in jlt_decode.py."""
import sys, os, time, functools, re, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import s70_field as F
import decode as D
import common as C
OUT = os.path.join(HERE, "records", "jlt"); os.makedirs(OUT, exist_ok=True)
F.photocurrent = functools.partial(F.photocurrent, resp=0.8, nep=15e-12)
K = 50; R = 0.01
mR = re.search(r'R(\d+)', sys.argv[1]) if len(sys.argv) > 1 else None
if mR:                                             # suffix R3, R5: reflectivity 3 %, 5 % (default 1 %)
    R = int(mR.group(1)) / 100.0
x_pm = np.linspace(-1240.0, 1240.0, 32)
CHIP_JLT = 511 / 16.352e-6
randi = np.random.default_rng(7).integers(0, 2, 511).astype(float)
mseq = F._mls(9)
z_uni = 2.5 * np.arange(1, K + 1)
G = np.exp(-(np.arange(1, 8) - 7.0) ** 2 / 8.0); G13 = np.r_[G, G[-2::-1]]; GT = np.r_[G13, G13]
det_case2 = np.zeros(K); det_case2[12:38] = 2 * 1.447 * 0.2e-9 * GT * 1e12
rng = np.random.default_rng(3)
for _ in range(1000):
    gaps = rng.uniform(2.2, 2.9, K - 1); gaps = gaps * (122.5 / np.sum(gaps))
    if gaps.min() >= 2.15:
        break
z_irr = 2.5 + np.r_[0.0, np.cumsum(gaps)]
band = np.repeat([0, 1, 2], [17, 17, 16])
det_bands = np.zeros(K)
for b in range(3):
    idx = np.flatnonzero(band == b); n = len(idx)
    det_bands[idx] = -800.0 + 800.0 * b + np.linspace(-350.0, 350.0, n)
def bands_gen(nb, pitch=44.0):
    """nb bands of K/nb gratings at a 44-pm pitch (as the three 0.8-nm bands), band centers spread over 2300/nb."""
    cnt = [K // nb + (1 if b < K % nb else 0) for b in range(nb)]
    d = np.zeros(K); i = 0
    for b in range(nb):
        c = (b - (nb - 1) / 2.0) * 2300.0 / nb
        d[i:i + cnt[b]] = c + pitch * (np.arange(cnt[b]) - (cnt[b] - 1) / 2.0); i += cnt[b]
    return d
jit = np.random.default_rng(11)
z_j5 = z_uni + jit.uniform(-0.05, 0.05, K)      # spacing tolerance +-5 cm, within the coherence length at 305 MHz
z_j50 = z_uni + jit.uniform(-0.5, 0.5, K)       # +-50 cm, beyond it
BASE = {
    'D0': (z_uni, np.zeros(K), CHIP_JLT, randi),
    'J5': (z_j5, np.zeros(K), CHIP_JLT, randi),
    'J50': (z_j50, np.zeros(K), CHIP_JLT, randi),
    'J5m': (z_j5, np.zeros(K), CHIP_JLT, mseq),       # tolerance +-5 cm, m-sequence
    'J5m50': (z_j5, np.zeros(K), 50e6, mseq),         # + 50 Mchip/s
    'D4j': (z_irr + np.random.default_rng(13).uniform(-0.05, 0.05, K), det_bands, 50e6, mseq),   # designed array with a +-5 cm position tolerance
    'J5mb': (z_j5, det_case2, CHIP_JLT, mseq),        # tolerance +-5 cm, m-sequence, case 2 of [JLT] (two bumps up to +0.58 nm)
    'D4x': (z_irr, det_bands + np.random.default_rng(5).uniform(-100.0, 100.0, K), 50e6, mseq),   # designed array, sensors detuned at random within +-100 pm
    'D4y': (z_irr, det_bands + np.random.default_rng(6).uniform(-200.0, 200.0, K), 50e6, mseq),   # detuned within the full +-200-pm range
    'D4g': (z_irr, det_bands + 100.0, 50e6, mseq),    # whole array shifted by +100 pm, the last grating at the end of the sweep
    'B2': (z_irr, bands_gen(2), 50e6, mseq),          # two bands of 25 gratings (R_c = 1/49)
    'B5': (z_irr, bands_gen(5), 50e6, mseq),          # five bands of 10 gratings
    'D1': (z_uni, np.zeros(K), CHIP_JLT, mseq),
    'D2': (z_uni, np.zeros(K), 50e6, mseq),
    'D3': (z_irr, np.zeros(K), 50e6, mseq),
    'D4': (z_irr, det_bands, 50e6, mseq),
}
name, seed = sys.argv[1], int(sys.argv[2])
mA = re.search(r'A(hn|hm|tk|ga|bl)', name)              # suffix Ahn/Ahm/Atk/Aga: grating apodization (default Blackman)
GT = (mA.group(1) if mA else 'bl') + '250'
name0 = re.sub(r'A(hn|hm|tk|ga|bl)', '', name)
name0 = re.sub(r'R\d+', '', name0)
direct = name0.endswith('d')                        # suffix d: direct paths only (ghosts=False)
name0 = name0[:-1] if direct else name0
core0 = re.sub(r'L\d+$', '', name0)
withrefs = core0.endswith('r')
core = core0.rstrip('r')
mL = re.search(r'L(\d+)$', name0)                     # suffix L1, L3, L10: source line 1, 3, 10 GHz (field oversampled accordingly)
wide = core.endswith('w')                           # legacy suffix w: 10-GHz line at the default oversampling
LW = float(mL.group(1)) * 1e9 if mL else (10e9 if wide else 305e6)
OSR = 4 if LW <= 1e9 else (8 if LW <= 3e9 else 32)
core = re.sub(r'L\d+$', '', core).rstrip('w')
z, det, chip, code = BASE[core]
z = np.sort(z)
sensors = [dict(z=float(z[k]), det=float(det[k]), R=R, g=GT) for k in range(K)]
refs = [dict(z=128.0, det=-950.0, R=R, g=GT), dict(z=131.0, det=950.0, R=R, g=GT)] if withrefs else []   # inside the sweep by more than a line width (references at +-1200 pm were cut by the sweep end and biased the correction)
out = OUT + "/jltc_%s_s%d.npz" % (name, seed)
t0 = time.time()
if not os.path.exists(out):
    F.run_array(sensors, refs, code=code, chip_rate=chip, x_pm=x_pm, seeds=(seed,), spc=64, linewidth_hz=LW, osr=OSR,
                laser_params=dict(eta=0.063), out=out, verbose=False, ghosts=not direct)
rec = np.load(out, allow_pickle=True)
spc = int(rec['spc']); ng = float(rec['n_group'])
bins = D.samples_of(rec['z'], ng, chip, spc)
X = np.array([D.correlate(rec['det'][0, m], code, spc) for m in range(len(x_pm))])
n = X.shape[1]
off = D.find_offset(X, bins[0], spc)
S = np.array([X[:, (int(round(b)) + off) % n] for b in bins]); scale = S[0].max(); S = S / scale
fit = np.array([C.gauss_fit_peak(x_pm, S[k], win_frac=1.0) for k in range(K)])
err = fit - det
line = '%-4s s%d %.0f s | raw: sigma %6.2f RMS %6.2f max %6.1f pm' % (name, seed, time.time() - t0, err.std(), np.sqrt(np.mean(err ** 2)), np.abs(err).max())
errc = err - err.mean()
line += ' | mean removed: sigma %6.2f' % errc.std()
if withrefs:
    bins_r = D.samples_of(rec['zr'], ng, chip, spc)
    Sr = np.array([X[:, (int(round(b)) + off) % n] for b in bins_r]) / scale
    fit_r = np.array([C.gauss_fit_peak(x_pm, Sr[j], win_frac=1.0) for j in range(len(bins_r))])
    err_r = fit_r - rec['detr_pm']
    p = np.polyfit(rec['detr_pm'], err_r, 1)
    errc = err - np.polyval(p, det)
    line += ' | 2 refs (offset+slope): sigma %6.2f RMS %6.2f max %6.1f pm' % (errc.std(), np.sqrt(np.mean(errc ** 2)), np.abs(errc).max())
line += ' | amp at 10,20,30,40,50: %s' % np.round(S.max(axis=1)[[9, 19, 29, 39, 49]], 2)
print(line, flush=True)
