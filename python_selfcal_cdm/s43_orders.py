"""s43_orders.py - how much the paths beyond the third order carry.

The ghost model of the paper stops at three reflections. A reviewer will ask
why, since the number of paths grows fast with the order: for an infinite
uniform array the count at each round trip is a Catalan number, split by the
number of reflections into Narayana numbers (Markowski et al. 2023, Sec. IV-A).

This enumerates every path of order 3 and 5 in a finite serial array of K
identical gratings, with the transmission of every grating the light passes
through, and compares the summed power of the two orders. A path of order n
visits gratings g1 > g2 < g3 > g4 < g5, in fiber order, and each crossing of a
grating costs (1 - R) in power.

The answer decides where the truncation is honest: at R = 1 percent it is,
at R = 10 percent it is only for short arrays, which is where the capacity
results already confine strong gratings.
"""
import numpy as np


def order3(K, R):
    a, b, c = np.meshgrid(np.arange(K), np.arange(K), np.arange(K),
                          indexing='ij')
    ok = (b < a) & (b < c)
    # crossings: out to a (a), back to b (a-b-1), out to c (c-b-1), home (c)
    cross = a + (a - b - 1) + (c - b - 1) + c
    p = R ** 3 * (1.0 - R) ** cross
    return int(ok.sum()), float(p[ok].sum())


def order5(K, R):
    idx = np.arange(K, dtype=np.int16)
    a, b, c, d, e = np.meshgrid(idx, idx, idx, idx, idx, indexing='ij')
    ok = (b < a) & (b < c) & (d < c) & (d < e)
    cross = (a + (a - b - 1) + (c - b - 1) + (c - d - 1) + (e - d - 1) + e)
    cross = cross.astype(np.float32)
    p = R ** 5 * (1.0 - R) ** cross
    return int(ok.sum()), float(p[ok].sum())


def direct(K, R):
    k = np.arange(K)
    return float((R * (1.0 - R) ** (2 * k)).sum())


print('%4s %6s %8s %10s %12s %12s %12s' % (
    'K', 'R', 'paths3', 'paths5', 'P5/P3', 'P3/Pdirect', 'P5/Pdirect'))
for K in (3, 8, 16, 32):
    for R in (0.10, 0.01):
        n3, p3 = order3(K, R)
        n5, p5 = order5(K, R)
        pd = direct(K, R)
        print('%4d %6.2f %8d %10d %12.4f %12.2e %12.2e'
              % (K, R, n3, n5, p5 / p3, p3 / pd, p5 / pd))
