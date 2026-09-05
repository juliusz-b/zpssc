"""s57_example_e2e.py - the eight-sensor temperature example of Table VII run
through the complete chain of s56, instead of the mechanism-by-mechanism
synthesis of s39.

Two arrays, same eight band positions (-175 .. +175 pm), sensors at 10 pm/K:
  initial   procured gratings, R = 10 percent, 250 pm, uniform 4 m, N = 127 at
            25 Mchip/s, two references, no deshadowing
  designed  R = 1 percent, 100 pm, positions 2.4 .. 33.1 m, N = 511 at
            100 Mchip/s, three references on a stub, deshadowing

Random per trial: line asymmetries, receiver noise, reference stability. The
chirp span is 0.2 FWHM in both. Reported per sensor: RMS over the trials, in
kelvin, plus the worst sensor and the count outside 1 K.

Output: out/s57_example_e2e.txt and .npz.
"""
import os
import numpy as np
import common as C
import s56_budget_e2e as S

PM_PER_K = C.TEMP_COEF_PM_PER_C
NU8 = np.linspace(-175.0, 175.0, 8)
NT = int(os.environ.get('S57_NTRIALS', '20'))

EX = {
    'initial':  dict(K=8, R=0.10, fwhm=250.0, nbits=7, B=25e6, spacing='uniform', d_min=4.0,
                     nref=2, peel=False, chirp=0.2, nub_pm=NU8),
    'designed': dict(K=8, R=0.01, fwhm=100.0, nbits=9, B=100e6, spacing='given',
                     z=[2.4, 5.1, 9.3, 12.2, 17.6, 21.0, 27.9, 33.1], d_min=2.0,
                     nref=3, peel=True, chirp=0.2, nub_pm=NU8),
}

if __name__ == '__main__':
    lines = []
    say = lambda s='': (print(s), lines.append(s))
    out = {}
    on = {k: True for k in S.ALL}
    for name, cfg in EX.items():
        for mode in ('baseline', 'offset'):
            E = np.array([S.one_trial(cfg, 2000 + t, on, mode) for t in range(NT)]) / PM_PER_K   # trials x 8, kelvin
            per = np.sqrt(np.mean(E ** 2, axis=0))
            say('=== %s array, %s receiver, %d trials' % (name, mode, NT))
            say('   per-sensor RMS [K]: ' + ' '.join('%.2f' % v for v in per))
            say('   worst sensor %.2f K, RMS over the eight %.2f K, sensors with RMS > 1 K: %d of 8, '
                'trial-wise worst |e| %.2f K' % (per.max(), np.sqrt(np.mean(per ** 2)), int((per > 1).sum()), np.abs(E).max()))
            out['%s__%s' % (name, mode)] = E
            # what dominates: without shadowing, without leakage
            for mech in ('shadow', 'leak', 'chirp', 'refs'):
                on_w = dict(on); on_w[mech] = False
                Ew = np.array([S.one_trial(cfg, 2000 + t, on_w, mode) for t in range(NT)]) / PM_PER_K
                say('   without %-7s RMS over the eight %.2f K' % (mech, np.sqrt(np.mean(Ew ** 2))))
            say()
    os.makedirs('out', exist_ok=True)
    open('out/s57_example_e2e.txt', 'w').write('\n'.join(lines) + '\n')
    np.savez('out/s57_example_e2e.npz', **out)
    print('saved out/s57_example_e2e.txt')
