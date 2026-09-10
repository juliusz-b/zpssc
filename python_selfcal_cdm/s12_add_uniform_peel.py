"""Dolicza do cache s12 (out/s12_results.npz) krzywa pojemnosci dla rownego odstepu z korekcja sekwencyjna
(cap_uni_peel), tymi samymi ziarnami co cap_uni. Uruchomic raz, potem S12_REPLOT=1 python s12_capacity.py."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')

src = open("s12_capacity.py", encoding="utf-8").read()
pre = src[:src.index("# (a) error vs array size at R = 10 percent, decomposed")]
pre = pre[:pre.rfind("# ----")]
g = {}
exec(compile(pre, "s12_prefix", "exec"), g)
sweep, TARGET_PM = g["sweep"], g["TARGET_PM"]
Rs = np.array([0.0003, 0.001, 0.003, 0.01, 0.03, 0.05, 0.10, 0.20, 0.30])
Ks_cap = np.array([4, 8, 16, 32, 48, 64, 96])


def capacity(R, mode, seed, **kw):
    e = sweep(Ks_cap, R, mode, 8, seed, **kw)
    if e[0] > TARGET_PM:
        return 0.0
    if e[-1] <= TARGET_PM:
        return float(Ks_cap[-1])
    i = int(np.argmax(e > TARGET_PM))
    x0, x1 = Ks_cap[i - 1], Ks_cap[i]
    y0, y1 = e[i - 1], e[i]
    return float(x0 + (TARGET_PM - y0) * (x1 - x0) / (y1 - y0))


cap = []
for i, R in enumerate(Rs):
    cap.append(capacity(R, 'uniform', 900 + i * 37, peel=True))
    print("R=%.4f  uniform+deshadowing K_max=%.1f" % (R, cap[-1]), flush=True)
cap_uni_peel = np.array(cap)
c = dict(np.load("out/s12_results.npz", allow_pickle=True))
c["cap_uni_peel"] = cap_uni_peel
np.savez("out/s12_results.npz", **c)
print("saved cap_uni_peel:", np.round(cap_uni_peel, 1))
