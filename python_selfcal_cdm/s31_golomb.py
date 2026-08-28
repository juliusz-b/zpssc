"""s31_golomb.py - grating placement as a Golomb ruler: zero ghost collisions.

Radar and array-processing aside (design session 2026-08-28): rule 1 of the
paper randomizes the grating pitch, which scatters about nine tenths of the
in-span third-order ghosts into empty delay bins. Antenna theory does the
same job deterministically. A Golomb ruler is a set of integer marks in which
every pairwise difference is distinct (minimum-redundancy arrays, sparse
rulers).

CLAIM. If the occupied delay bins form a Golomb ruler, NO third-order ghost
lands on an occupied bin. Proof in one line: a collision means
tau_a - tau_b + tau_c = tau_d, i.e. tau_a - tau_b = tau_d - tau_c, two equal
pairwise differences, which the ruler forbids unless a = d and b = c, and
b = c is not a valid path (the intermediate grating must differ from the
last reflector).

Bonus for this paper: the optimal 14-mark Golomb ruler has length exactly
127, so an N = 127 code hosts a perfect 14-grating layout with the ruler
filling the delay span completely.

Cost: a K-mark ruler needs about K^2 bins, so the trick covers K up to
about sqrt(N-ish) gratings per band (14 at N = 127, 25 at N = 511 does NOT
fit: L(25) = 480 <= 511 does fit, checked below if marks known). For larger
K randomized spacing remains the tool.

Output: fraction of third-order ghost paths landing on occupied bins for
uniform, randomized, and Golomb layouts, plus the count of in-span ghosts.
"""
import numpy as np

# known optimal rulers (Shearer's tables), verified for the Golomb property
RULERS = {
    10: [0, 1, 6, 10, 23, 26, 34, 41, 53, 55],
    14: [0, 4, 6, 20, 35, 52, 59, 77, 78, 86, 89, 99, 122, 127],
}


def is_golomb(marks):
    m = np.asarray(marks)
    d = (m[None, :] - m[:, None])[np.triu_indices(len(m), 1)]
    return len(np.unique(d)) == len(d)


def ghost_stats(bins, span):
    """All third-order paths (a, b, c): b strictly before a and c, c != b.
    Returns (paths, in-span ghosts, collisions with occupied bins)."""
    s = set(int(x) for x in bins)
    bins = sorted(s)
    paths = inspan = coll = 0
    for ib, b in enumerate(bins):
        later = bins[ib + 1:]
        for a in later:
            for c in later:
                paths += 1
                g = a - b + c
                if g <= span:
                    inspan += 1
                    if g in s:
                        coll += 1
    return paths, inspan, coll


N = 127
K = 14
span = N

print('=== third-order ghost collisions, K = %d gratings, %d delay bins ===' % (K, N))

# uniform pitch, same span
u = np.arange(K) * (N // (K - 1))
p, i, c = ghost_stats(u, span)
print('  uniform pitch      : %4d paths, %4d in span, %4d on occupied bins (%.0f%%)'
      % (p, i, c, 100.0 * c / max(i, 1)))

# randomized, 200 trials
rng = np.random.default_rng(31)
fr = []
for _ in range(200):
    r = np.sort(rng.choice(np.arange(1, N + 1), size=K, replace=False))
    p, i, c = ghost_stats(r, span)
    fr.append(c / max(i, 1))
print('  randomized (200 tr): mean %.1f%% of in-span ghosts on occupied bins,'
      ' best trial %.1f%%' % (100 * np.mean(fr), 100 * np.min(fr)))

# Golomb ruler
g14 = RULERS[14]
assert is_golomb(g14), 'ruler 14 fails the Golomb check'
assert max(g14) == 127, 'ruler 14 should span exactly 127'
p, i, c = ghost_stats(g14, span)
print('  golomb ruler K=14  : %4d paths, %4d in span, %4d on occupied bins (%.0f%%)'
      % (p, i, c, 100.0 * c / max(i, 1)))
print('  (optimal 14-mark ruler length = 127 = N: the code span is used whole)')

print()
print('=== the same at K = 10 (ruler length 55, loose in 127 bins) ===')
g10 = [x + 1 for x in RULERS[10]]
assert is_golomb(g10)
p, i, c = ghost_stats(g10, span)
print('  golomb ruler K=10  : %4d paths, %4d in span, %4d collisions' % (p, i, c))
u10 = np.arange(10) * 14 + 1
p, i, c = ghost_stats(u10, span)
print('  uniform K=10       : %4d paths, %4d in span, %4d collisions (%.0f%%)'
      % (p, i, c, 100.0 * c / max(i, 1)))

print()
print('VERDICT: a Golomb-ruler layout removes ghost-bin collisions exactly,')
print('not statistically, for K up to the ruler that fits the code span.')
print('Randomized spacing stays the fallback for larger K.')
