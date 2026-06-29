"""S5: odpornosc kotwicy referencyjnej na wielokrotne odbicia (ghosty), ghost
przesuniety widmowo wzgledem primary. Dopasowanie Gaussa do dominujacego primary
jest mocno odporne na pociaganie przez ghost; przy identycznych siatkach pozycja
jest scisle niezmiennicza (skaluje sie tylko amplituda)."""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore'); import common as C
fwhm=15.0; nu=np.linspace(-120,160,3001); PM=C.PM_PER_GHZ
def peak(alpha,s,asym=0.0):
    prim=C.fbg_tanh(nu,0.0,fwhm,asym)
    ghost=C.fbg_tanh(nu,s,fwhm,asym) if s>0 else C.fbg_tanh(nu,0.0,fwhm,asym)
    S=prim+alpha*ghost
    p_gauss=C.gauss_fit_peak(nu,S)
    p_cent=C.centroid(nu,S,thr=0.05)   # naiwny centroid calego widma
    return p_gauss*PM, p_cent*PM, S.max()
alphas=np.linspace(0,0.20,11)  # realistyczny zakres ghost/primary (po korelacji)
g_id=[peak(a,0.0)[0] for a in alphas]
g_s=[peak(a,40.0)[0] for a in alphas]        # ghost +0.32 nm
c_s=[peak(a,40.0)[1] for a in alphas]        # naiwny centroid (dla kontrastu)
amp=[peak(a,40.0)[2] for a in alphas]
Rs=np.linspace(0.02,0.30,15); ratio=Rs**2   # ghost/primary ~ R^2 po korelacji
gR=[abs(peak(rr,40.0)[0]) for rr in ratio]       # dopasowanie Gaussa (uzywane)
cR=[abs(peak(rr,40.0)[1]) for rr in ratio]       # naiwny centroid (dla kontrastu)
fig,ax=plt.subplots(1,2,figsize=(10.4,3.6))
ax[0].plot(alphas,g_id,'o-',label='Gauss, siatki identyczne')
ax[0].plot(alphas,g_s,'s-',label='Gauss, ghost +0.32 nm')
ax[0].plot(alphas,c_s,'^--',color='0.5',label='naiwny centroid (ghost +0.32 nm)')
ax[0].set_xlabel('amplituda ghost / primary'); ax[0].set_ylabel('przesuniecie pozycji [pm]')
ax[0].set_xlim(0,0.2); ax[0].set_title('(a) Pozycja kotwicy: Gauss odporny vs centroid'); ax[0].legend(fontsize=7,loc='upper left')
ax[1].plot(Rs*100,cR,'^--',color='0.5',label='naiwny centroid')
ax[1].plot(Rs*100,gR,'o-',label='dopasowanie Gaussa (uzywane)')
ax[1].set_xlabel('reflektancja siatki R [%]'); ax[1].set_ylabel('blad kotwicy z ghostow [pm]')
ax[1].set_title('(b) Blad kotwicy vs R (ghost ratio ~ R$^2$)'); ax[1].legend(fontsize=7,loc='upper left')
plt.tight_layout(); plt.savefig('figs/fig_s5_ghosts.png',dpi=140)
print(f"Gauss @alpha=0.2: identyczne={g_id[-1]:.2f} pm, ghost+0.32nm={g_s[-1]:.2f} pm; naiwny centroid={c_s[-1]:.1f} pm")
print(f"blad kotwicy Gauss @R=9%: {abs(peak(0.09**2,40.0)[0]):.3f} pm")
print('saved figs/fig_s5_ghosts.png')
