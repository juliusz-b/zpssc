"""S3 (GLOWNY WYNIK): rezyduum bledu dlugosci fali po korekcie kanalami
referencyjnymi, w funkcji stosunku chirp/FWHM i niedopasowania nachylenia
(ksztaltu) miedzy siatka referencyjna a czujnikowa.
Chirp dzieli sie na czesc wspolna (offset zmienny w pasmie -> usuwalny przez
referencje) i czesc znieksztalcenia ksztaltu (zalezna od nachylenia siatki ->
rezyduum, ktore zostaje po odjeciu trybu wspolnego)."""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore'); import common as C
np.random.seed(3)
fwhm=15.0; band=600.0   # GHz (~10 nm)
def mean_off(nu, delta):  # adiabatyczny offset chirpu: rosnie z Δ i zmienny w pasmie (tryb wspolny)
    return delta*(0.20 + 0.20*(nu/band) + 0.15*(nu/band)**2)
def app_err(nub, delta, asym, skew=1.2):
    nu=np.linspace(nub-5*fwhm, nub+5*fwhm, 1401)
    spec=C.fmam_readout(nu, nub, fwhm, delta, mean_off=mean_off(nub,delta), skew=skew, shape='tanh', asym=asym)
    return (C.gauss_fit_peak(nu, spec)-nub)*C.PM_PER_GHZ   # pm
def calibrate(se, sx, re, rx, order):
    if order<0: return se
    order=min(order, len(rx)-1)
    p=np.polyfit(rx, re, order)
    return se - np.polyval(p, sx)
nref=4; nsen=12
ref_nu=np.linspace(-band*0.8, band*0.8, nref)
sen_nu=np.sort(np.random.uniform(-band*0.85, band*0.85, nsen))
ref_asym=np.random.uniform(-0.05,0.05,nref)   # referencje: dobrze scharakteryzowane
sen_asym=np.random.uniform(-0.30,0.30,nsen)   # czujniki: rozne nachylenia/ksztalty
deltas=np.linspace(0.05,1.0,12)*fwhm
res={}
for label,nr,order in [("brak kalibracji",0,-1),("1 ref (offset)",1,0),
                        ("2 ref (liniowo)",2,1),("4 ref (kwadratowo)",4,2)]:
    rr=[]
    for d in deltas:
        se=np.array([app_err(sen_nu[i],d,sen_asym[i]) for i in range(nsen)])
        if nr>0:
            re=np.array([app_err(ref_nu[j],d,ref_asym[j]) for j in range(nr)])
            cor=calibrate(se,sen_nu,re,ref_nu[:nr],order)
        else:
            cor=se
        rr.append(np.sqrt(np.mean(cor**2)))
    res[label]=rr
# panel (a): blad w pasmie przy Δ=0.5·FWHM, surowy vs po 4-ref
d_fix=0.5*fwhm
se=np.array([app_err(sen_nu[i],d_fix,sen_asym[i]) for i in range(nsen)])
re=np.array([app_err(ref_nu[j],d_fix,ref_asym[j]) for j in range(nref)])
cor=calibrate(se,sen_nu,re,ref_nu,2)
# panel (c): rezyduum vs |niesymetria czujnika| (proxy niedopasowania nachylenia)
d_c=0.6*fwhm; asyms=np.linspace(0,0.4,9); rva=[]
for a in asyms:
    se2=np.array([app_err(sen_nu[i],d_c,a*(1 if sen_asym[i]>=0 else -1)) for i in range(nsen)])
    re2=np.array([app_err(ref_nu[j],d_c,ref_asym[j]) for j in range(nref)])
    rva.append(np.sqrt(np.mean(calibrate(se2,sen_nu,re2,ref_nu,2)**2)))
fig,ax=plt.subplots(1,3,figsize=(13.2,3.6))
ax[0].axhline(0,c='0.8',lw=0.7)
ax[0].plot(sen_nu/1000*C.PM_PER_GHZ, se,'o-',label='surowy blad')
ax[0].plot(sen_nu/1000*C.PM_PER_GHZ, cor,'s-',label='po korekcie (4 ref)')
ax[0].plot(ref_nu/1000*C.PM_PER_GHZ, re,'k^',ms=8,label='referencje')
ax[0].set_xlabel('pozycja w pasmie [nm]'); ax[0].set_ylabel(r'blad $\lambda_B$ [pm]')
ax[0].set_title('(a) Blad w pasmie, Δ=0.5·FWHM'); ax[0].legend(fontsize=7)
for label,rr in res.items(): ax[1].plot(deltas/fwhm, rr,'o-',ms=4,label=label)
ax[1].set_xlabel('chirp Δ / FWHM siatki'); ax[1].set_ylabel('rezyduum RMS [pm]')
ax[1].set_title('(b) Rezyduum vs chirp i liczba referencji'); ax[1].legend(fontsize=7)
ax[2].plot(asyms, rva,'o-')
ax[2].set_xlabel('niesymetria / nachylenie czujnika'); ax[2].set_ylabel('rezyduum RMS [pm]')
ax[2].set_title('(c) Rezyduum vs niedopasowanie nachylenia')
plt.tight_layout(); plt.savefig('figs/fig_s3_residual.png',dpi=140)
print(f"surowy RMS @Δ=0.5FWHM = {np.sqrt(np.mean(se**2)):.1f} pm;  po 4-ref = {np.sqrt(np.mean(cor**2)):.1f} pm")
for label,rr in res.items():
    print(f"  {label:18s}: @Δ=0.2·FWHM={rr[2]:.1f} pm | @Δ=0.5·FWHM={rr[5]:.1f} pm | @Δ=FWHM={rr[-1]:.1f} pm")
print("saved figs/fig_s3_residual.png")
