"""Przelicza cztery krzywe pojemnosci z panelu (b) s12 z wieksza liczba ukladow (24 zamiast 8) i zapisuje je
do cache out/s12_results.npz (klucze cap_uni, cap_rnd, cap_uni_peel, cap_peel, cap_trials).
Potem: S12_REPLOT=1 python s12_capacity.py."""
import sys
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')

NT = int(sys.argv[1]) if len(sys.argv) > 1 else 24
src = open("s12_capacity.py", encoding="utf-8").read()
pre = src[:src.index("# (a) error vs array size at R = 10 percent, decomposed")]
pre = pre[:pre.rfind("# ----")]
g = {}
exec(compile(pre, "s12_prefix", "exec"), g)
sweep, TARGET_PM = g["sweep"], g["TARGET_PM"]
Rs = np.array([0.0003, 0.001, 0.003, 0.01, 0.03, 0.05, 0.10, 0.20, 0.30])
Ks_cap = np.array([4, 8, 16, 32, 48, 64, 96])


def capacity(R, mode, seed, **kw):
    e = sweep(Ks_cap, R, mode, NT, seed, **kw)
    if e[0] > TARGET_PM:
        return 0.0
    if e[-1] <= TARGET_PM:
        return float(Ks_cap[-1])
    i = int(np.argmax(e > TARGET_PM))
    x0, x1 = Ks_cap[i - 1], Ks_cap[i]
    y0, y1 = e[i - 1], e[i]
    return float(x0 + (TARGET_PM - y0) * (x1 - x0) / (y1 - y0))


t0 = time.time()
out = {}
for key, mode, base, kw in (("cap_uni", "uniform", 900, {}), ("cap_rnd", "random", 500, {}),
                            ("cap_uni_peel", "uniform", 900, {"peel": True}), ("cap_peel", "random", 500, {"peel": True})):
    out[key] = np.array([capacity(R, mode, base + i * 37, **kw) for i, R in enumerate(Rs)])
    print("%s done %.0f s: %s" % (key, time.time() - t0, np.round(out[key], 1)), flush=True)
c = dict(np.load("out/s12_results.npz", allow_pickle=True))
c.update(out)
c["cap_trials"] = np.array(NT)
np.savez("out/s12_results.npz", **c)
print("saved, trials =", NT, flush=True)
