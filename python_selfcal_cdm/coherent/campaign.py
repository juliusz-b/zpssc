"""campaign.py - records of the coherent model for the bound check of Fig. 10(b).

python campaign.py capR [R]     K = 4, 8, 16, 32 gratings at R = 1, 3, 10, 30 %, two random layouts each
                                (minimum spacing 4 m), 300-MHz line, all reflection orders
python campaign.py capRb [R]    the same at a 1-GHz line and 8-m minimum spacing (filled markers of Fig. 10b)
python campaign.py capRd        the R = 1 % and 3 % layouts of capR without multiple reflections (ghosts=False),
                                to isolate the coherent ghost term
Optional second argument: one reflectivity (0.01, 0.03, 0.1, 0.3). Records go to records/, layouts to
records/*_layout.npz. One record of 32 gratings takes about a minute. fig_bound.py decodes them.
"""
import os
import sys
import time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import runs as VV

X64 = VV.X64
LW = 300e6
what = sys.argv[1]
SEL = float(sys.argv[2]) if len(sys.argv) > 2 else None
t0 = time.time()


def layout(rng, K, dmin_m):
    """K positions with a minimum spacing dmin_m in random bins (as in the power model, s12/s13)."""
    return 2.0 + np.cumsum(dmin_m * (1.0 + rng.uniform(0.0, 1.0, K)))


def tick(msg):
    print("%s %.0f s" % (msg, time.time() - t0), flush=True)


if what == "capR":
    rng = np.random.default_rng(91)
    for R in (0.01, 0.03, 0.10, 0.30):
        if SEL is not None and abs(R - SEL) > 1e-6:
            continue
        for K in (4, 8, 16, 32):
            for t in range(2):
                zz = layout(rng, K, 4.0); dd = rng.uniform(-200.0, 200.0, K)
                np.savez(os.path.join(VV.OUT, "capR_R%02d_K%03d_t%d_layout.npz" % (int(R * 100), K, t)), z=zz, det=dd)
                VV.run("capR_R%02d_K%03d_t%d" % (int(R * 100), K, t), VV.sens(zz, dd, R, "bl250"), x_pm=X64, linewidth_hz=LW)
            tick("capR R=%.2f K=%d" % (R, K))
elif what == "capRb":
    rng = np.random.default_rng(91)
    for R in (0.01, 0.03, 0.10, 0.30):
        if SEL is not None and abs(R - SEL) > 1e-6:
            continue
        for K in (4, 8, 16, 32):
            for t in range(2):
                zz = layout(rng, K, 8.0); dd = rng.uniform(-200.0, 200.0, K)
                VV.run("capRb_R%02d_K%03d_t%d" % (int(R * 100), K, t), VV.sens(zz, dd, R, "bl250"), x_pm=X64, linewidth_hz=1000e6)
            tick("capRb R=%.2f K=%d" % (R, K))
elif what == "capRd":
    for R, K in ((0.01, 8), (0.01, 16), (0.03, 8)):
        for t in (0, 1):
            lay = np.load(os.path.join(VV.OUT, "capR_R%02d_K%03d_t%d_layout.npz" % (int(R * 100), K, t)))
            VV.run("capRd_R%02d_K%03d_t%d" % (int(R * 100), K, t), VV.sens(lay["z"], lay["det"], R, "bl250"), x_pm=X64, linewidth_hz=LW, ghosts=False)
            tick("capRd R=%.2f K=%d t=%d" % (R, K, t))
else:
    raise SystemExit("unknown campaign: " + what)
print("done", what, "%.0f s" % (time.time() - t0))
