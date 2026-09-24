"""runs.py - helpers to run s70_field.run_array and cache the records in records/, plus the Rule A check.

python runs.py ruleA     two gratings (R = 10 %, 4 and 16 m), shift of grating 2 against the detuning of grating 1,
                         300-MHz line. Writes cache/ruleA_shift.npz (Fig. 3a via fig_validation.py) and prints the
                         table against the closed-form shadowing shift (eq:lawA).
Records already present in records/ are reused, so a run can be resumed.
"""
import os
import sys
import time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import common as C
import s70_field as F
import decode as D

OUT = os.path.join(HERE, "records")      # field sampled 4 times denser than the detector (osr=4)
os.makedirs(OUT, exist_ok=True)
X64 = np.linspace(-650.0, 650.0, 64)
NU8 = np.linspace(-175.0, 175.0, 8)
C0 = F.C0


def sens(zs, dets, R, g):
    return [dict(z=float(z), det=float(d), R=R, g=g) for z, d in zip(zs, dets)]


def refs_of(n, R, g, stub=(60.0, 75.0, 90.0)):
    dets = np.linspace(-180.0, 180.0, n) if n > 1 else [0.0]
    return [dict(z=stub[j], det=float(dets[j]), R=R, g=g) for j in range(n)]


def run(name, sensors, refs=(), **kw):
    p = os.path.join(OUT, name + ".npz")
    if os.path.exists(p):
        return np.load(p, allow_pickle=True)
    t0 = time.time()
    F.run_array(sensors, refs, out=p, **kw)
    print("   %s: %.0f s" % (name, time.time() - t0), flush=True)
    return np.load(p, allow_pickle=True)


def width_of(z, k, seed=0):
    x = z["x_pm"]; code = z["code"]; spc = int(z["spc"]); ng = float(z["n_group"]); rate = float(z["chip_rate"])
    X = np.array([D.correlate(z["det"][seed, m], code, spc) for m in range(len(x))])
    off = D.find_offset(X, D.samples_of(z["z"][k], ng, rate, spc), spc)
    return C.gauss_fit_full(x, D.read_nearest(X, D.samples_of(z["z"][k], ng, rate, spc), off, spc))[2]


def centre_of(z, k, mode="baseline", seed=0):
    x = z["x_pm"]; code = z["code"]; spc = int(z["spc"]); ng = float(z["n_group"]); rate = float(z["chip_rate"])
    X = np.array([D.correlate(z["det"][seed, m], code, spc, mode) for m in range(len(x))])
    off = D.find_offset(X, D.samples_of(z["z"][k], ng, rate, spc), spc)      # filter delay from the read grating
    return C.gauss_fit_peak(x, D.read_nearest(X, D.samples_of(z["z"][k], ng, rate, spc), off, spc))


def v_ruleA(lw=300e6):
    """Rule A: two gratings of 10 % at 4 and 16 m, error of grating 2 against the detuning of grating 1."""
    tag = "_lw%03d" % int(lw / 1e6)
    ref = run("ruleA_D3000" + tag, sens([4.0, 16.0], [3000.0, 0.0], 0.10, "bl250"), x_pm=X64, linewidth_hz=lw)
    c0 = centre_of(ref, 1)
    sig = width_of(ref, 1)
    print("fitted sigma of the line in the coherent model: %.0f pm (grating FWHM 250 pm)" % sig)
    print("| Delta [pm] | coherent model [pm] | shadowing shift (eq:lawA) [pm] |")
    print("|---|---|---|")
    ds = np.array([-260, -130, 0, 50, 100, 130, 170, 200, 260, 330, 400])
    err = []
    for Dl in ds:
        z = run("ruleA_D%+04d%s" % (Dl, tag), sens([4.0, 16.0], [float(Dl), 0.0], 0.10, "bl250"), x_pm=X64, linewidth_hz=lw)
        c = centre_of(z, 1) - c0
        f = -4.0 / 3.0 * np.sqrt(2.0 / 3.0) * 0.10 * Dl * np.exp(-Dl ** 2 / (3 * sig ** 2))
        print("| %+d | %+.1f | %+.1f |" % (Dl, c, f))
        err.append(c)
    np.savez(os.path.join(HERE, "cache", "ruleA_shift.npz"), detuning_pm=ds, shift_pm=np.array(err), sigma_pm=sig, R=0.10, linewidth_hz=lw)
    print("saved cache/ruleA_shift.npz")


if __name__ == "__main__":
    names = sys.argv[1:] or ["ruleA"]
    for n in names:
        globals()["v_" + n]()
