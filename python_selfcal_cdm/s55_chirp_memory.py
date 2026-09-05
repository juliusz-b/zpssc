"""s55_chirp_memory.py - chirp that depends on the bit history, and the echoes it
makes in the correlation.

A chip that follows a zero carries the transient overshoot of the laser, a chip
that follows a one does not. The received waveform is therefore not the code
times one constant. For an m-sequence the chips that follow a zero,
v_n = u_n (1 - u_{n-1}), are a sum of shifted copies of the same sequence
(shift-and-add), so their correlation with the replica has peaks of +-1/2 at
lag 0, lag 1 and one further lag set by the generator polynomial. A grating
that sits one chip (or that further lag) closer to the source than the grating
being read adds an echo shaped like the derivative of its line, and the fitted
centre moves by about half the mean overshoot within a chip. This does not
shrink with N and the code autocorrelation does not show it.

Time-domain model: unipolar m-sequence, SPC samples per chip, xi(t) = adiabatic
offset while on + overshoot A_T exp(-t/tau) after every 0->1 edge. Two gratings,
k under test at delay 0 and a neighbour j at DD chips (negative = closer to the
source), detuned by DL. Mean removal, periodic correlation with the bipolar
replica, chip-average sampling, 64-step sweep, Gaussian LS fit of the centre.
The reference case replaces the transient by its power-weighted mean on every
on chip, which is what a chirp-kernel model does. The difference is the memory
effect.

Output: out/s55_chirp_memory.txt (no figure).
"""
import os, warnings
import numpy as np
from scipy.signal import max_len_seq
import common as C
warnings.filterwarnings('ignore')

SPC = 32
FWHM = 250.0
SIG = FWHM / 2.35482
LAM = np.linspace(-650, 650, 64)


class Model:
    def __init__(self, nbits):
        self.N = 2 ** nbits - 1
        self.u = max_len_seq(nbits)[0].astype(float)
        self.b = 2 * self.u - 1
        self.drive = np.repeat(self.u, SPC)
        self.L = len(self.drive)

    def echo_lags(self):
        v = self.u * (1 - np.roll(self.u, 1))
        beta = np.array([(np.roll(self.b, l) * v).sum() for l in range(self.N)]) / (self.u * self.b).sum()
        return [(int(l), round(float(beta[l]), 2)) for l in np.where(np.abs(beta) > 0.1)[0]]

    def xi(self, A_ad, A_T, tau_chip, memory=True):
        x = np.full(self.L, A_ad, float)
        for e in np.where(np.diff(np.r_[self.u[-1], self.u]) > 0)[0]:
            i0 = e * SPC
            x[i0:] += A_T * np.exp(-np.arange(self.L - i0) / (tau_chip * SPC))
        if memory:
            return x
        return np.full(self.L, (x * self.drive).sum() / self.drive.sum())

    def spectrum(self, xi, DD, DL):
        sh = int(round(DD * SPC))
        xi_j, dr_j = np.roll(xi, sh), np.roll(self.drive, sh)
        out = np.zeros(len(LAM))
        for m, lam in enumerate(LAM):
            y = self.drive * np.exp(-0.5 * ((lam + xi) / SIG) ** 2) \
                + dr_j * np.exp(-0.5 * ((lam - DL + xi_j) / SIG) ** 2)
            y -= y.mean()
            out[m] = (y.reshape(self.N, SPC).mean(axis=1) * self.b).sum()
        return out

    def shift(self, A_ad, A_T, tau_chip, DD, DL):
        c1 = C.gauss_fit_peak(LAM, self.spectrum(self.xi(A_ad, A_T, tau_chip, True), DD, DL))
        c0 = C.gauss_fit_peak(LAM, self.spectrum(self.xi(A_ad, A_T, tau_chip, False), DD, DL))
        return float(c1 - c0)

    def worst(self, A_ad, A_T, tau_chip, DD):
        return max((self.shift(A_ad, A_T, tau_chip, DD, DL) for DL in (0.0, 100.0, 200.0)), key=abs)


lines = []
def say(s=''):
    print(s); lines.append(s)

say('pattern-echo lags of the scipy m-sequences (lag, beta):')
for nb in (6, 7, 9, 10):
    say('  N=%4d: %s' % (2 ** nb - 1, Model(nb).echo_lags()))

M = Model(7)
say('')
say('N=127, FWHM=250 pm, adiabatic offset 20 pm. Centre shift memory - no memory [pm],')
say('worst over neighbour detuning 0/100/200 pm. DD<0: neighbour closer to the source.')
DDs = (-1, -7, -2, -3, -13, 1, 7)
say('%-26s' % 'overshoot A_T, tau [chip]' + ''.join('%8s' % ('DD=%d' % d) for d in DDs) + '   alone')
for A_T, tau in ((30.0, 0.1), (30.0, 0.3), (30.0, 1.0), (80.0, 0.3), (80.0, 1.0)):
    row = [M.worst(20.0, A_T, tau, d) for d in DDs]
    alone = M.shift(20.0, A_T, tau, 50, 1e9)
    say('%-26s' % ('%3.0f pm, %.1f' % (A_T, tau)) + ''.join('%8.2f' % x for x in row) + '%8.2f' % alone)

say('')
say('bench: three gratings 4 m apart, overshoot 30 pm with tau = 12 ns')
for B in (25e6, 50e6, 100e6):
    chip_m = 3e8 / 1.468 / B / 2           # one chip of round-trip delay in metres
    DD = -4.0 / chip_m
    tau = 12e-9 * B
    if abs(DD - round(DD)) < 1e-6:
        DD = int(round(DD))
    w = max((M.shift(20.0, 30.0, tau, DD, DL) for DL in (0.0, 100.0, 200.0)), key=abs)
    say('  %3.0f Mchip/s: 4 m = %.1f chips, tau = %.2f chip, neighbour one pitch closer: %6.2f pm'
        % (B / 1e6, -DD, tau, w))

os.makedirs('out', exist_ok=True)
open('out/s55_chirp_memory.txt', 'w').write('\n'.join(lines) + '\n')
print('saved out/s55_chirp_memory.txt')
