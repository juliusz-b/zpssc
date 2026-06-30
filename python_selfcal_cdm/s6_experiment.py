"""s6_experiment.py - bench-realistic prediction for the procured hardware.
Three identical FBGS DTG-A3A4 gratings (order S-2026-0066), all at ~1545 nm,
FWHM ~250 pm, R=10%, in series and separated only by code/delay (CDM). Each
grating sits on a Peltier stage (~10 pm/C) so its Bragg wavelength is tunable
over ~+/-225 pm. Two gratings act as co-coded references at fixed set-points; the
third is the sensor, swept across the tuning band. The script predicts (i) that
the three spectrally-overlapping gratings separate by code, and (ii) the apparent
Bragg error from source-chirp FM-AM and how far co-coded references suppress it,
at the scale actually achievable on this bench (3 gratings used, 15 procured)."""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore'); import common as C
np.random.seed(6)
F  = C.FBG_FWHM_GHZ                      # ~31.25 GHz (= 250 pm)
PM = C.PM_PER_GHZ
R  = C.FBG_R                            # 0.10
# Draw-tower gratings are very uniform; small per-grating lineshape differences:
asym = np.array([+0.04, -0.03, +0.02])

def app_err(nu_b, delta, a, skew=1.2):
    """Apparent Bragg error [pm] of one grating under source-chirp FM-AM."""
    g = np.linspace(nu_b - 6*F, nu_b + 6*F, 1400)
    spec = C.fmam_readout(g, nu_b, F, delta, mean_off=0.30*delta, skew=skew, shape='tanh', asym=a)
    return (C.gauss_fit_peak(g, spec) - nu_b) * PM

# (a) three overlapping spectra at three temperature set-points
dnu = np.array([C.temp_to_dnu_ghz(dT) for dT in (-12.0, 0.0, +18.0)])   # Celsius -> GHz
nu = np.linspace(-80, 80, 2000)
specs = [C.fbg_tanh(nu, dnu[k], F, asym[k]) for k in range(3)]
composite = np.clip(sum(R*s for s in specs), 0, None)

# (b) sensor apparent error vs chirp excursion, raw vs 2-reference corrected
ref_nu = np.array([-18.0, +18.0])       # two references set ~+/-150 pm by temperature
ref_as = np.array([+0.01, -0.01])       # references well characterized
sen_nu = 4.0; sen_as = 0.05             # sensor near band center
deltas = np.linspace(0.3, 8.0, 16)      # chirp excursion [GHz] (set by modulation depth)
raw=[]; cor=[]
for d in deltas:
    e_s = app_err(sen_nu, d, sen_as)
    e_r = np.array([app_err(ref_nu[j], d, ref_as[j]) for j in range(2)])
    p = np.polyfit(ref_nu, e_r, 1)      # linear correction from 2 references
    raw.append(abs(e_s)); cor.append(abs(e_s - np.polyval(p, sen_nu)))

# (c) temperature calibration at a fixed realistic modulation (chirp 4 GHz)
d_fix = 4.0
true_nu = np.linspace(-25, 25, 21)      # sensor true Bragg via temperature [GHz]
e_r = np.array([app_err(ref_nu[j], d_fix, ref_as[j]) for j in range(2)])
p = np.polyfit(ref_nu, e_r, 1)
e_raw = np.array([app_err(tn, d_fix, sen_as) for tn in true_nu])
e_cor = np.array([app_err(tn, d_fix, sen_as) - np.polyval(p, tn) for tn in true_nu])

# (d) CDM despreading of the 3 overlapping gratings (the bench measurement)
M=128; codes=C.gold_set(7); N=codes.shape[1]
slots=np.linspace(-70,70,M); taus=[7,29,61]      # distinct delays (4 m spacing)
r=np.zeros(N)
for k in range(3):
    refl=C.fbg_tanh(slots,dnu[k],F,asym[k])*R
    for mi,a_m in enumerate(refl):
        r+=a_m*np.roll(codes[mi%N],taus[k])
desp=[np.array([C.periodic_xcorr(r,codes[mi%N])[taus[k]] for mi in range(M)]) for k in range(3)]
dmax=max(d.max() for d in desp)
comp_slots=np.clip(sum(C.fbg_tanh(slots,dnu[k],F,asym[k])*R for k in range(3)),0,None)

fig,ax=plt.subplots(2,2,figsize=(12,7.2))
ax[0,0].plot(nu*PM/1000, composite,'k',lw=1.6,label='composite (overlapping)')
for k in range(3): ax[0,0].plot(nu*PM/1000, R*specs[k],'--',lw=1,label=f'FBG{k+1} (T set-point)')
ax[0,0].set_xlabel('wavelength offset from 1545 nm [nm]'); ax[0,0].set_ylabel('reflectivity')
ax[0,0].set_title('(a) 3 identical FBGs @1545 nm, FWHM 250 pm, R=10%'); ax[0,0].legend(fontsize=7)
ax[0,1].plot(deltas/F, raw,'o-',label='raw (no calibration)')
ax[0,1].plot(deltas/F, cor,'s-',label='2 co-coded references')
ax[0,1].axvspan(0.3/F, 4.0/F, color='0.9'); ax[0,1].text(2.2/F,max(raw)*0.5,'realistic\nmodulation',fontsize=7,ha='center')
ax[0,1].set_xlabel('chirp excursion / FWHM'); ax[0,1].set_ylabel('|apparent Bragg error| [pm]')
ax[0,1].set_title('(b) Sensor error vs chirp (FWHM = 250 pm)'); ax[0,1].legend(fontsize=7)
ax[1,0].axhline(0,c='0.85',lw=0.7)
ax[1,0].plot(true_nu*PM/1000, e_raw,'o-',label='raw')
ax[1,0].plot(true_nu*PM/1000, e_cor,'s-',label='reference-corrected')
ax[1,0].plot(ref_nu*PM/1000,[0,0],'k^',ms=10,label='references (fixed T)')
ax[1,0].set_xlabel('sensor true Bragg offset (temperature) [nm]'); ax[1,0].set_ylabel('apparent error [pm]')
ax[1,0].set_title('(c) Temperature calibration, chirp = 4 GHz'); ax[1,0].legend(fontsize=7)
ax[1,1].plot(slots*PM/1000, comp_slots/comp_slots.max(),'k',lw=1,label='composite (input)')
for k in range(3): ax[1,1].plot(slots*PM/1000, desp[k]/dmax,'.-',ms=3,label=f'despread FBG{k+1}')
ax[1,1].set_xlabel('wavelength offset [nm]'); ax[1,1].set_ylabel('norm. amplitude')
ax[1,1].set_title('(d) CDM despreading separates the 3 overlapping FBGs'); ax[1,1].legend(fontsize=7)
plt.tight_layout(); plt.savefig('figs/fig_s6_experiment.png',dpi=140)
i4=int(np.argmin(abs(deltas-4)))
print(f"FBG: lambda0={C.LAMBDA0_NM} nm, FWHM={C.FBG_FWHM_PM} pm (={F:.1f} GHz), R={C.FBG_R}, {C.N_GRATINGS_BENCH} on bench / {C.N_GRATINGS_AVAIL} procured")
print(f"raw sensor error @chirp=4 GHz (Delta/FWHM={4/F:.2f}): {raw[i4]:.1f} pm; 2-ref corrected: {cor[i4]:.2f} pm")
print(f"temperature-calibration residual RMS (corrected): {np.sqrt(np.mean(e_cor**2)):.2f} pm")
print("saved figs/fig_s6_experiment.png")
