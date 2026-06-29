"""S1: wlasciwosci korelacyjne kodow Gold/Kasami i podloga listkow bocznych.
Listki boczne (deterministyczne) wyznaczaja podloge rozplotu CDM."""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import common as C
g=C.gold_set(7)
c0=g[0]; auto=C.periodic_xcorr(c0,c0); cross=C.periodic_xcorr(g[0],g[3])
ga,gx=C.code_corr_stats(g)
print(f"Gold N=127: max auto-sidelobe={ga:.3f}, max cross={gx:.3f}, teoria 17/127={17/127:.3f}")
fig,ax=plt.subplots(1,2,figsize=(9.2,3.3))
L=len(c0); lag=np.arange(L)-L//2
ax[0].plot(lag,np.roll(auto,L//2),lw=0.9,label='autokorelacja')
ax[0].plot(lag,np.roll(cross,L//2),lw=0.9,label='korelacja wzajemna')
ax[0].axhline(17/127,ls='--',c='k',lw=0.7); ax[0].axhline(-17/127,ls='--',c='k',lw=0.7,label='granica Golda 17/127')
ax[0].set_xlabel('opoznienie [chip]'); ax[0].set_ylabel('korelacja (znorm.)')
ax[0].set_title('(a) Kody Gold, N=127'); ax[0].legend(fontsize=7)
ns=[5,6,7]; fl=[]
for n in ns:
    gg=C.gold_set(n); a,x=C.code_corr_stats(gg); fl.append(max(a,x))
ax[1].semilogy(ns,fl,'o-',label='podloga (max listek/cross)')
ax[1].semilogy(ns,[1/np.sqrt(2**n-1) for n in ns],'s--',label=r'$1/\sqrt{N}$')
ax[1].set_xlabel('stopien kodu n'); ax[1].set_ylabel('podloga korelacji')
ax[1].set_title('(b) Podloga rozplotu vs dlugosc kodu'); ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig('figs/fig_s1_codes.png',dpi=140)
print("saved figs/fig_s1_codes.png")
