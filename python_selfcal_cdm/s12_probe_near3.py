"""Sonda wokol R = 3 %: pojemnosc z korekcja dla rownego i losowego odstepu przy R = 2..4 %, te same odstrojenia
(jedno ziarno) dla wszystkich R, K = 32, 48, 64, 24 ukladow."""
import sys
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')

src = open("s12_capacity.py", encoding="utf-8").read()
pre = src[:src.index("# (a) error vs array size at R = 10 percent, decomposed")]
pre = pre[:pre.rfind("# ----")]
g = {}
exec(compile(pre, "s12_prefix", "exec"), g)
sweep, TARGET_PM = g["sweep"], g["TARGET_PM"]
Ks = np.array([32, 48, 64])
NT = 24
t0 = time.time()
print("K =", Ks, " trials =", NT, flush=True)
for R in (0.02, 0.025, 0.03, 0.035, 0.04):
    for mode, seed in (("uniform", 1048), ("random", 648)):
        e = sweep(Ks, R, mode, NT, seed, peel=True)
        i = int(np.argmax(e > TARGET_PM)) if (e > TARGET_PM).any() else None
        if i is None:
            cap = "> %d" % Ks[-1]
        elif i == 0:
            cap = "< 32"
        else:
            cap = "%.1f" % (Ks[i-1] + (TARGET_PM - e[i-1]) * (Ks[i] - Ks[i-1]) / (e[i] - e[i-1]))
        print("R=%.3f %-7s +desh  RMS %s  K_10pm %s   (%.0f s)" % (R, mode, np.round(e, 1), cap, time.time() - t0), flush=True)
print("done", flush=True)
