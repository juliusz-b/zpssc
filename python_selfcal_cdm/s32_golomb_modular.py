"""s32_golomb_modular.py - the Golomb-ruler claim under a PERIODIC code.

Referee finding (29.08.2026): the s31 argument ignored two things. First,
the optimal 14-mark ruler has LENGTH 127, so it occupies bins 0..127, which
is 128 bins, one more than an N = 127 code offers, and the last mark aliases
onto the first. Second, the correlation is periodic, so a ghost delay beyond
the span does not vanish, it folds back modulo N and can land on an occupied
bin. Both points are checked here, and the repair is checked too: a MODULAR
Golomb ruler (all pairwise differences distinct modulo N) makes the collision
argument exact under wrap-around. Its size is bounded near sqrt(N).

Also counted: near-collisions, ghosts landing exactly one bin from an
occupied mark, inside the two-chip triangular correlation footprint.
"""
import numpy as np

RULER_A = [0, 4, 6, 20, 35, 52, 59, 77, 78, 86, 89, 99, 122, 127]   # s31
RULER_B = [0, 5, 28, 38, 41, 49, 50, 68, 75, 92, 107, 121, 123, 127]  # referee


def ghosts(marks):
    """All third-order paths (a, b, c) with b strictly before a and c,
    c != b (a == c allowed). Returns list of raw ghost delays."""
    m = sorted(set(int(x) for x in marks))
    out = []
    for ib, b in enumerate(m):
        later = m[ib + 1:]
        for a in later:
            for c in later:
                out.append(a - b + c)
    return np.array(out)


def stats(marks, N):
    occ = set(int(x) % N for x in marks)
    g = ghosts(marks)
    gw = g % N
    beyond = int(np.sum(g >= N))
    coll = int(np.sum([x in occ for x in gw]))
    near = int(np.sum([((x + 1) % N in occ) or ((x - 1) % N in occ) for x in gw]))
    return len(g), beyond, coll, near


def is_modular_golomb(marks, N):
    m = np.asarray(marks)
    d = (m[None, :] - m[:, None]) % N
    d = d[~np.eye(len(m), dtype=bool)]
    return len(np.unique(d)) == len(d)


def max_marks(N):
    """Counting bound on a modular Golomb ruler (Sidon set) mod N.

    Every ordered pair of distinct marks gives a nonzero difference mod N,
    and all k(k-1) of them must differ, so k(k-1) <= N-1. This is a theorem,
    not a search result, and it is what the paper should quote.
    """
    k = 1
    while (k + 1) * k <= N - 1:
        k += 1
    return k


def search_modular(N, k, tries, rng):
    """Randomized greedy search for a k-mark modular Golomb ruler mod N.

    Greedy with many restarts beats depth-first search here: the tree is
    wide, a node budget truncates it long before it pays off, and a fresh
    random order explores more of the space per unit time.
    """
    for _ in range(tries):
        marks = [0]
        used = set()
        for x in rng.permutation(np.arange(1, N)):
            d = {(x - y) % N for y in marks} | {(y - x) % N for y in marks}
            if len(d) == 2 * len(marks) and not (d & used):
                marks.append(int(x))
                used |= d
                if len(marks) == k:
                    return sorted(marks)
    return None


rng = np.random.default_rng(32)
N = 127

print('=== 14-mark optimal rulers under a periodic N = 127 code ===')
for name, r in (('ruler A (s31)', RULER_A), ('ruler B (referee)', RULER_B)):
    print('  %-18s distinct bins mod N: %d of 14  (mark 127 == mark 0)'
          % (name, len(set(x % N for x in r))))
    p, beyond, coll, near = stats(r, N)
    print('      paths %d, beyond span %d, collisions after wrap %d (%.1f%%),'
          ' one bin from a mark %d' % (p, beyond, coll, 100.0 * coll / p, near))

print()
print('=== references under the same wrap-around rule, K = 14 ===')
u = np.arange(14) * 9
p, beyond, coll, near = stats(u, N)
print('  uniform pitch 9  : collisions %d of %d (%.0f%%)' % (coll, p, 100.0 * coll / p))
fr = []
for _ in range(300):
    r = np.sort(rng.choice(np.arange(N), size=14, replace=False))
    p, beyond, coll, near = stats(r, N)
    fr.append(coll / p)
print('  randomized (300) : mean %.1f%% collisions' % (100 * np.mean(fr)))

print()
print('=== repair: modular Golomb rulers (all differences distinct mod N) ===')
for NN in (63, 127, 255, 511, 1023):
    kb = max_marks(NN)
    best = None
    for k in range(6, kb + 1):
        m = search_modular(NN, k, tries=3000, rng=rng)
        if m is None:
            break
        best = (k, m)
    k, m = best
    assert is_modular_golomb(m, NN)
    p, beyond, coll, near = stats(m, NN)
    print('  N = %4d: counting bound k <= %2d, search found %2d, collisions after wrap %d of %d,'
          ' within one bin %d' % (NN, kb, k, coll, p, near))
    if NN == 127:
        print('            marks at N = 127:', m)
