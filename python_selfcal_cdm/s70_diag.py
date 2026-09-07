"""Diagnostics of the rate-equation laser of s70_field: chirp, extinction, record boundary."""
import time
import numpy as np
import s70_field as F

fs = 25e6 * 64
I = F.drive_current(F._mls(7), 64, 0.033, 0.018)
t0 = time.time()
E, f_c = F.laser_record(I, fs, seed=1)
print("laser %.1f s, carrier offset f_c = %.3f GHz" % (time.time() - t0, f_c / 1e9))
ph = np.unwrap(np.angle(E)); fi = np.r_[np.diff(ph), 0] / (2 * np.pi) * fs
P = np.abs(E) ** 2; on = I > 0.033; off = I < 0.033; pm = F.LAM0 ** 2 / F.C0 * 1e12
print("residual inst freq [GHz]: on %.3f, off %.3f, on-off %.3f = %.1f pm, std within on %.3f" % (fi[on].mean() / 1e9, fi[off].mean() / 1e9, (fi[on].mean() - fi[off].mean()) / 1e9, (fi[on].mean() - fi[off].mean()) * pm, fi[on].std() / 1e9))
print("power on/off %.2f; boundary |E|^2/mean first 3 %s last 3 %s" % (P[on].mean() / P[off].mean(), np.round(P[:3] / P.mean(), 2), np.round(P[-3:] / P.mean(), 2)))
chip = P.reshape(127, 64)
print("mean chip power profile (16 bins) of on-chips after zero:", np.round(chip[np.array(F._mls(7)) * (1 - np.roll(F._mls(7), 1)) > 0].mean(axis=0).reshape(16, -1).mean(axis=1) / P.max(), 2))
