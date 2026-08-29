"""s33_ncpc_robustness.py - does the +1 replica survive fractional delays,
jitter and drift? (referee finding B3, 29.08.2026)

The referee claims that a return displaced by a fraction f of a chip leaks
f/(1-f) of its peak into EVERY other bin under the +1 replica, so the scheme
loses to the -1/N floor already at f = 0.008 chip. That would be fatal,
because rule 3 randomizes grating positions. It is checked here with an
oversampled waveform: rectangular chips, fractional delay, boxcar integration
over one chip, then chip-rate samples. Linearity says the answer in advance:
a fractional shift is a two-tap mix of adjacent integer shifts, and the
cyclic inverse maps each to a delta, so the peak splits between two adjacent
bins and nothing reaches the others. Random per-sample jitter and a linear
intensity ramp are checked the same way, for both receivers.
"""
import numpy as np
import common as C

N = 127
OS = 200                                    # samples per chip, oversampled
c01 = C._mls01(7).astype(float)
msb = 2.0 * c01 - 1.0                       # +1 replica (sum = +1)
wave = np.repeat(c01, OS)                   # one period, rectangular chips


def record(delay_chips, ramp=0.0, jitter=0.0, rng=None):
    """Chip-rate samples of one period returned with a (fractional) delay,
    boxcar-integrated over each chip, optional linear ramp of relative
    amplitude 'ramp' across the period and per-sample timing jitter."""
    shift = int(round(delay_chips * OS))
    w = np.roll(wave, shift)
    if ramp:
        w = w * (1.0 + ramp * (np.arange(N * OS) / (N * OS) - 0.5))
    if jitter and rng is not None:
        jit = rng.normal(0.0, jitter * OS, N).astype(int)
        idx = (np.arange(N) * OS + OS // 2 + jit) % (N * OS)   # sample at chip centre
        return w[idx]                        # instantaneous samples with jitter
    return w.reshape(N, OS).mean(axis=1)


def rx_plus(y):
    return np.real(np.fft.ifft(np.fft.fft(y) * np.fft.fft(msb).conj()))


def rx_mean(y):
    return np.real(np.fft.ifft(np.fft.fft(y - y.mean()) * np.fft.fft(msb).conj()))


def leak(z, k):
    """Largest |output| outside bins k and k+1, relative to the peak."""
    pk = np.abs(z[[k, (k + 1) % N]]).max()
    mask = np.ones(N, bool); mask[[k, (k + 1) % N]] = False
    return np.abs(z[mask]).max() / pk


print('=== fractional delay f (chips), boxcar sampling, no noise ===')
print('   f       +1 replica leak   mean-removed leak   (floor -1/N = %.4f)' % (1.0 / N))
for f in (0.0, 0.05, 0.1, 0.25, 0.5):
    y = record(40 + f)
    print('  %4.2f     %10.2e         %10.4f' % (f, leak(rx_plus(y), 40), leak(rx_mean(y), 40)))

print()
print('=== random per-sample timing jitter (RMS, chips), 50 trials ===')
rng = np.random.default_rng(33)
for sj in (0.02, 0.05, 0.1, 0.2):
    lp, lm = [], []
    for _ in range(50):
        y = record(40, jitter=sj, rng=rng)
        lp.append(leak(rx_plus(y), 40)); lm.append(leak(rx_mean(y), 40))
    print('  sigma_t = %.2f: +1 replica %.4f, mean-removed %.4f' % (sj, np.mean(lp), np.mean(lm)))

print()
print('=== linear intensity ramp across one period ===')
for r in (0.0, 0.05, 0.127, 0.3):
    y = record(40, ramp=r)
    print('  ramp %.3f: +1 replica %.4f, mean-removed %.4f' % (r, leak(rx_plus(y), 40), leak(rx_mean(y), 40)))

print()
print('=== constant offset d (per unit peak) ===')
for d in (0.1, 0.5):
    y = record(40) + d
    zp = rx_plus(y); z0 = rx_plus(record(40))
    ped = zp - z0
    print('  d = %.1f: pedestal mean %.4f, std %.1e (flat -> baseline term removes it)'
          % (d, ped.mean() / z0[40], ped.std() / z0[40]))
