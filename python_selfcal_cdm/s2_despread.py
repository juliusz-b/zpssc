"""S2: CDM despreading of spectrally overlapping gratings + reference channels.
Each of M wavelength slots carries its own code; gratings are separated by delay.
Shows that (i) CDM separates overlapping spectra and (ii) reference channels are
simply additional CDM channels (the same despreading path)."""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')
import common as C
np.random.seed(1)
codes=C.gold_set(7)             # 129 codes, N=127
N=codes.shape[1]
M=128                           # wavelength slots (dense: ~4 slots / FWHM)
band=(-260.0,260.0)
lam=np.linspace(band[0],band[1],M)
fwhm=15.0
gratings=[("S1",-150.0,5),("R1",-70.0,17),("S2",20.0,33),("R2",110.0,55),("S3",190.0,83)]
R=0.09
r=np.zeros(N)
for (name,nub,tau) in gratings:
    Rk=R*C.fbg_gauss(lam,nub,fwhm)
    for m in range(M):
        r += Rk[m]*np.roll(codes[m],tau)
r += 0.0015*np.random.randn(N)
composite=sum(R*C.fbg_gauss(lam,nub,fwhm) for (_,nub,_) in gratings)
corr=np.array([C.periodic_xcorr(r,codes[m]) for m in range(M)])   # [M,N]
fig,ax=plt.subplots(1,2,figsize=(9.6,3.4))
ax[0].plot(lam,composite,'k',lw=1.5); ax[0].set_title('(a) Composite spectrum (uncoded system)')
ax[0].set_xlabel('detuning [GHz]'); ax[0].set_ylabel('reflection (sum)')
errs=[]
for (name,nub,tau) in gratings:
    spec=corr[:,tau]; est=C.gauss_fit_peak(lam,spec); errs.append((est-nub)*C.PM_PER_GHZ)
    ls='-' if name.startswith('S') else '--'
    ax[1].plot(lam,spec,ls,lw=1.2,label=name)
ax[1].set_title('(b) After CDM despreading (separated channels)'); ax[1].set_xlabel('detuning [GHz]')
ax[1].set_ylabel('recovered reflection'); ax[1].legend(fontsize=7,ncol=2)
plt.tight_layout(); plt.savefig('figs/fig_s2_despread.png',dpi=140)
print("lambda read-out errors [pm]:",{g[0]:round(e,1) for g,e in zip(gratings,errs)})
print(f"MAE (ideal channel, no chirp) = {np.mean(np.abs(errs)):.1f} pm")
af,cf=C.code_corr_stats(codes)
print(f"sidelobe/cross floor = {max(af,cf):.3f} (Gold N=127)")
print("saved figs/fig_s2_despread.png")
