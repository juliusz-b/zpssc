"""s7_axis_nonlinearity.py - the nonlinear voltage->wavelength tuning of the
HCG-VCSEL, how co-coded reference gratings straighten that axis, and the role of
the Peltier budget.

Physics. The interrogator ramps the tuning voltage v over [0,1]. The TRUE optical
frequency is nonlinear in v (HCG-VCSEL ~ quadratic, BW10-1550):
    nu_true(v) = common.vcsel_sweep(v, span, quad, cube).
A naive interrogator assumes a straight line nu_lin(v) = v*span. A grating with
true Bragg nu_B reflects at the voltage v* where nu_true(v*) = nu_B; the naive
readout reports nu_lin(v*), so the axis error is nu_lin(v*) - nu_B. Reference
gratings at temperature-set, KNOWN frequencies give anchors (v*, nu_B): 2 refs ->
linear axis fit (leaves curvature), 3 refs -> quadratic fit, MZI k-clock -> dense
true map. On top sits the source-chirp FM-AM (common.fmam_readout at the measured
chirp); the MZI is blind to it, so FM-AM always needs the references.

Peltier budget. Only REFERENCES need a Peltier (known, stabilized wavelength);
SENSORS do not. With common.N_PELTIER = 3 stages the array carries at most 3
references (plus unlimited non-Peltier sensors). This script contrasts two
regimes: (top) a wide WDM-CDM band where the axis nonlinearity is large and 3
references or an MZI are needed; (bottom) the procured narrow bench (all gratings
within ~+/-225 pm of 1545 nm by temperature) where the local nonlinearity is mild
and 2-3 references suffice, so 3 Peltiers are enough.
"""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore'); import common as C

PM=C.PM_PER_GHZ; F=C.FBG_FWHM_GHZ
SPAN=600.0; QUAD=0.30; CUBE=0.10; CHIRP=C.MEAS_CHIRP_GHZ; M=2048
v_grid=np.linspace(0.0,1.0,M)
nu_true=C.vcsel_sweep(v_grid,SPAN,quad=QUAD,cube=CUBE)
nu_lin=v_grid*SPAN
v_star=lambda nb: float(np.interp(nb,nu_true,v_grid))
nu_naive=lambda nb: float(np.interp(v_star(nb),v_grid,nu_lin))
def fmam_err(nb,asym,skew=1.2):
    g=np.linspace(nb-6*F,nb+6*F,1200)
    spec=C.fmam_readout(g,nb,F,CHIRP,mean_off=0.30*CHIRP,skew=skew,shape='tanh',asym=asym)
    return C.gauss_fit_peak(g,spec)-nb
mae=lambda e: np.mean(np.abs(e))*PM
LABELS=['none','2 ref','2 ref\n+ MZI','3 ref','MZI only','3 ref\n+ MZI']

def run_band(nu_lo, nu_hi, ref_fracs, nsen=12, nr=30, seed=7):
    """MAE per strategy for gratings in the optical-frequency window [nu_lo,nu_hi]."""
    ref_nu=nu_lo+np.array(ref_fracs)*(nu_hi-nu_lo)
    rng=np.random.default_rng(seed); acc=np.zeros(6); last=None
    for it in range(nr):
        sen_nu=np.sort(rng.uniform(nu_lo,nu_hi,nsen))
        sen_as=rng.uniform(-0.30,0.30,nsen); ref_as=rng.uniform(-0.05,0.05,len(ref_nu))
        sv=np.array([v_star(n) for n in sen_nu]); rv=np.array([v_star(n) for n in ref_nu])
        ax_sen=np.array([nu_naive(n) for n in sen_nu])-sen_nu
        fm_sen=np.array([fmam_err(sen_nu[i],sen_as[i]) for i in range(nsen)])
        fm_ref=np.array([fmam_err(ref_nu[j],ref_as[j]) for j in range(len(ref_nu))])
        # axis residual after K-reference polyfit of nu_B vs v*
        axres=lambda K,o: np.polyval(np.polyfit(rv[:K],ref_nu[:K],o),sv)-sen_nu
        fmres=lambda K,o: fm_sen-np.polyval(np.polyfit(ref_nu[:K],fm_ref[:K],o),sen_nu)
        e_none = ax_sen+fm_sen
        e_2    = axres(2,1)+fmres(2,1)            # 2 ref: linear axis + linear FM-AM
        e_2mzi = fmres(2,1)                        # 2 ref + MZI: axis exact, FM-AM via 2 ref
        e_3    = axres(3,2)+fmres(3,2)            # 3 ref: quadratic axis + quadratic FM-AM
        e_mzi  = fm_sen                            # MZI only: axis exact, FM-AM raw
        e_3mzi = fmres(3,2)                        # 3 ref + MZI: axis exact, FM-AM via 3 ref
        acc+=np.array([mae(e_none),mae(e_2),mae(e_2mzi),mae(e_3),mae(e_mzi),mae(e_3mzi)])
        if it==0: last=(sen_nu,ref_nu,e_none,e_2,e_3,e_3mzi)
    return acc/nr, last

# wide WDM-CDM band (gratings written across ~4.8 nm) and narrow procured bench
wide_vals,wide_last = run_band(0.05*SPAN, 0.95*SPAN, [0.05,0.50,0.95])
BENCH_W=56.0                                   # +/-225 pm window (= +/-28 GHz) ...
bc=0.16*SPAN                                   # ... placed near the 1545 nm edge of the sweep
bench_vals,bench_last = run_band(bc-BENCH_W/2, bc+BENCH_W/2, [0.10,0.50,0.90])

print(f"HCG sweep span={SPAN:.0f} GHz, quad={QUAD}, cube={CUBE}; chirp={CHIRP:.1f} GHz; N_PELTIER={C.N_PELTIER}")
print("WIDE WDM-CDM band (gratings across ~4.8 nm):")
for l,v in zip(LABELS,wide_vals): print(f"  {l.replace(chr(10),' '):12s}: {v:7.1f} pm")
print(f"NARROW procured bench (+/-225 pm window):")
for l,v in zip(LABELS,bench_vals): print(f"  {l.replace(chr(10),' '):12s}: {v:7.1f} pm")

# ---- figure: top = wide band, bottom = procured bench ----
fig,ax=plt.subplots(2,2,figsize=(12,8))
cols=['#c0392b','#e67e22','#16a085','#2980b9','#8e44ad','#27ae60']
def pos_panel(a,last,title):
    sn,rn,e0,e2,e3,e3m=last
    a.axhline(0,c='0.85',lw=0.7)
    a.plot(sn*PM/1000,e0*PM,'o-',ms=3,color=cols[0],label='none')
    a.plot(sn*PM/1000,e2*PM,'s-',ms=3,color=cols[1],label='2 ref (linear)')
    a.plot(sn*PM/1000,e3*PM,'D-',ms=3,color=cols[3],label='3 ref (quadratic)')
    a.plot(sn*PM/1000,e3m*PM,'^-',ms=3,color=cols[5],label='3 ref + MZI')
    a.plot(rn*PM/1000,np.zeros(len(rn)),'k^',ms=9,label='references (Peltier)')
    a.set_xlabel('position in band [nm]'); a.set_ylabel(r'apparent $\lambda_B$ error [pm]')
    a.set_title(title); a.legend(fontsize=6.5,ncol=2)
def bar_panel(a,vals,title):
    a.bar(range(6),vals,color=cols)
    a.set_xticks(range(6)); a.set_xticklabels(LABELS,fontsize=7)
    a.set_ylabel('MAE [pm] (mean of 30 runs)'); a.set_title(title)
    for i,v in enumerate(vals): a.text(i,v+max(vals)*0.02,f'{v:.0f}',ha='center',fontsize=8)
pos_panel(ax[0,0],wide_last,'(a) Wide WDM-CDM band: apparent error vs position')
bar_panel(ax[0,1],wide_vals,'(b) Wide band: MAE per strategy (axis nonlinearity dominates)')
pos_panel(ax[1,0],bench_last,'(c) Procured bench (+/-225 pm): apparent error vs position')
bar_panel(ax[1,1],bench_vals,'(d) Bench: MAE per strategy (3 Peltiers suffice)')
plt.tight_layout(); plt.savefig('figs/fig_s7_axis.png',dpi=140); print("saved figs/fig_s7_axis.png")
