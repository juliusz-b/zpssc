"""S4: comparison of calibration strategies (mandatory baseline), averaged over
realizations. Wavelength axis: sweep nonlinearity (smooth + higher-order wiggle)
+ drift, plus chirp FM-AM. The MZI k-clock straightens the AXIS (nonlinearity +
drift) but does NOT remove the grating-edge FM-AM; the reference channels remove
FM-AM in common mode, but a few points cannot reproduce the nonlinearity wiggle.
Only the combination removes both."""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore'); import common as C
fwhm=15.0; band=600.0; delta=0.5*fwhm; PM=C.PM_PER_GHZ
NLA=50.0; WIG=8.0; drift=2.0
def nl(nu): return NLA*((nu/band)-(nu/band)**3)+WIG*np.sin(3*np.pi*nu/band)
def mean_off(nu,d): return d*(0.20+0.20*(nu/band)+0.15*(nu/band)**2)
def fmam_err(nub,asym,skew=1.2):
    nu=np.linspace(nub-5*fwhm,nub+5*fwhm,1001)
    spec=C.fmam_readout(nu,nub,fwhm,delta,mean_off=mean_off(nub,delta),skew=skew,shape='tanh',asym=asym)
    return (C.gauss_fit_peak(nu,spec)-nub)
nref=4; nsen=12; NR=40
ref_nu=np.linspace(-band*0.8,band*0.8,nref)
mae=lambda e: np.mean(np.abs(e))*PM
acc=np.zeros(4); last=None
rng=np.random.default_rng(4)
for it in range(NR):
    ref_asym=rng.uniform(-0.05,0.05,nref)
    sen_nu=np.sort(rng.uniform(-band*0.85,band*0.85,nsen)); sen_asym=rng.uniform(-0.3,0.3,nsen)
    ax_ref=nl(ref_nu)+drift; ax_sen=nl(sen_nu)+drift
    fm_ref=np.array([fmam_err(ref_nu[j],ref_asym[j]) for j in range(nref)])
    fm_sen=np.array([fmam_err(sen_nu[i],sen_asym[i]) for i in range(nsen)])
    tot_ref=ax_ref+fm_ref; tot_sen=ax_sen+fm_sen
    e_nocal=tot_sen; e_mzi=fm_sen                       # MZI removes axis, FM-AM remains
    e_ref=tot_sen-np.polyval(np.polyfit(ref_nu,tot_ref,2),sen_nu)
    e_refmzi=fm_sen-np.polyval(np.polyfit(ref_nu,fm_ref,2),sen_nu)
    acc+=np.array([mae(e_nocal),mae(e_mzi),mae(e_ref),mae(e_refmzi)])
    if it==0: last=(sen_nu,e_nocal,e_mzi,e_ref,e_refmzi)
vals=acc/NR
labels=['no\ncalibration','MZI k-clock\nonly','references\nonly','references\n+ MZI']
print({l.replace(chr(10),' '):round(v,1) for l,v in zip(labels,vals)})
fig,ax=plt.subplots(1,2,figsize=(10.6,3.8))
ax[0].bar(range(4),vals,color=['#c0392b','#e67e22','#2980b9','#27ae60'])
ax[0].set_xticks(range(4)); ax[0].set_xticklabels(labels,fontsize=8)
ax[0].set_ylabel(f'MAE [pm] (mean of {NR} runs)'); ax[0].set_title('(a) MAE vs calibration strategy')
for i,v in enumerate(vals): ax[0].text(i,v+1.5,f'{v:.0f}',ha='center',fontsize=9)
sn,e0,e1,e2,e3=last
ax[1].axhline(0,c='0.85',lw=0.7)
ax[1].plot(sn/1000*PM,e0*PM,'o-',ms=4,label='none')
ax[1].plot(sn/1000*PM,e1*PM,'^-',ms=4,label='MZI')
ax[1].plot(sn/1000*PM,e2*PM,'s-',ms=4,label='references')
ax[1].plot(sn/1000*PM,e3*PM,'D-',ms=4,label='ref+MZI')
ax[1].set_xlabel('position in band [nm]'); ax[1].set_ylabel('error [pm]')
ax[1].set_title('(b) Per-sensor error (one run)'); ax[1].legend(fontsize=7,ncol=2)
plt.tight_layout(); plt.savefig('figs/fig_s4_baselines.png',dpi=140); print('saved figs/fig_s4_baselines.png')
