"""Exhaustive search for a modular Golomb ruler (cyclic Sidon set) with K marks modulo N.

All K(K-1) ordered differences must be distinct and nonzero mod N. The counting bound is
K(K-1) <= N-1. The search fixes the marks 0 and 1, which loses nothing whenever some
difference of the set is a unit mod N (always true for prime N): multiply the set by the
inverse of that difference and translate. Marks are enumerated in ascending order.

Result used in the paper (Section on ghost-free spacing):
    N=127, K=11  ->  no ruler exists, the largest has 10 marks (134e6 nodes, ~20 min).
Sanity cases: N=13 K=4 -> {0,1,3,9}, N=31 K=6, N=73 K=9, N=91 K=10 all found (Singer sets).

Usage: python s54_modular_ruler.py N K [time_limit_s]
"""
import sys
import time


def search(N, K, limit_s=600.0):
    t0 = time.time()
    used = [False] * N
    marks = [0, 1]
    used[1] = used[N - 1] = True
    best = [2]
    nodes = [0]

    def rec():
        nodes[0] += 1
        if (nodes[0] & 0xFFFFF) == 0 and time.time() - t0 > limit_s:
            raise TimeoutError
        depth = len(marks)
        best[0] = max(best[0], depth)
        if depth == K:
            return True
        for x in range(marks[-1] + 1, N - (K - depth) + 1):
            ok = True
            diffs = []
            for m in marks:
                d1 = (x - m) % N
                d2 = N - d1
                if used[d1] or used[d2] or d1 == d2:
                    ok = False
                    break
                used[d1] = used[d2] = True
                diffs.append(d1)
                diffs.append(d2)
            if ok:
                marks.append(x)
                if rec():
                    return True
                marks.pop()
            for d in diffs:
                used[d] = False
        return False

    try:
        found = rec()
    except TimeoutError:
        return None, best[0], nodes[0], time.time() - t0
    return (list(marks) if found else False), best[0], nodes[0], time.time() - t0


if __name__ == "__main__":
    N, K = int(sys.argv[1]), int(sys.argv[2])
    limit = float(sys.argv[3]) if len(sys.argv) > 3 else 600.0
    res, best, nodes, dt = search(N, K, limit)
    if res is None:
        status = "time limit reached"
    elif res is False:
        status = "no ruler"
    else:
        status = "found %s" % res
    print("N=%d K=%d: %s (largest depth %d, %d nodes, %.0f s)" % (N, K, status, best, nodes, dt))
