"""S5: immunity of the reference anchor to multiple reflections (ghosts); the
ghost is spectrally offset from the primary. The Gaussian fit to the dominant
primary is strongly resistant to ghost pull; for identical gratings the position
is exactly invariant (only the amplitude scales)."""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore'); import common as C
fwhm=15.0; nu=np.linspace(-120,160,3001); PM=C.PM_PER_GHZ
def peak(alpha,s,asym=0.0):
    prim=C.fbg_tanh(nu,0.0,fwhm,asym)
    ghost=C.fbg_tanh(nu,s,fwhm,asym) if s>0 else C.fbg_tanh(nu,0.0,fwhm,asym)
    S=prim+alpha*ghost
    return C.gauss_fit_peak(nu,S)*PM, C.centroid(nu,S,thr=0.05)*PM, S.max()
alphas=np.linspace(0,0.20,11)  # realistic ghost/primary range (after correlation)
g_id=[peak(a,0.0)[0] for a in alphas]
g_s=[peak(a,40.0)[0] for a in alphas]        # ghost +0.32 nm
c_s=[peak(a,40.0)[1] for a in alphas]        # naive centroid (for contrast)
Rs=np.linspace(0.02,0.30,15); ratio=Rs**2   # ghost/primary ~ R^2 after correlation
gR=[abs(peak(rr,40.0)[0]) for rr in ratio]       # Gaussian fit (used)
cR=[abs(peak(rr,40.0)[1]) for rr in ratio]       # naive centroid (for contrast)
fig,ax=plt.subplots(1,2,figsize=(10.4,3.6))
ax[0].plot(alphas,g_id,'o-',label='Gauss, identical gratings')
ax[0].plot(alphas,g_s,'s-',label='Gauss, ghost +0.32 nm')
ax[0].plot(alphas,c_s,'^--',color='0.5',label='naive centroid (ghost +0.32 nm)')
ax[0].set_xlabel('ghost / primary amplitude'); ax[0].set_ylabel('position shift [pm]')
ax[0].set_xlim(0,0.2); ax[0].set_title('(a) Anchor position: Gauss robust vs centroid'); ax[0].legend(fontsize=7,loc='upper left')
ax[1].plot(Rs*100,cR,'^--',color='0.5',label='naive centroid')
ax[1].plot(Rs*100,gR,'o-',label='Gaussian fit (used)')
ax[1].set_xlabel('grating reflectivity R [%]'); ax[1].set_ylabel('anchor error from ghosts [pm]')
ax[1].set_title('(b) Anchor error vs R (ghost ratio ~ R$^2$)'); ax[1].legend(fontsize=7,loc='upper left')
plt.tight_layout(); plt.savefig('figs/fig_s5_ghosts.png',dpi=140)
print(f"Gauss @alpha=0.2: identical={g_id[-1]:.2f} pm, ghost+0.32nm={g_s[-1]:.2f} pm; naive centroid={c_s[-1]:.1f} pm")
print(f"Gauss anchor error @R=9%: {abs(peak(0.09**2,40.0)[0]):.3f} pm")
print('saved figs/fig_s5_ghosts.png')
