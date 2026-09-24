"""Fig. 3 and Fig. 11 of the paper from the coherent model.

Fig. 3 (fig_validation_curvature): (a) direct return of the 8th co-tuned grating from coupled-mode spectra,
    crossing the critical reflectivity R_c = 1/15, against the Gaussian-line approximation (dotted),
    (b) the same eight gratings of unequal width in the two orders, (c) all 40320 orders (cache/order_perm.npz,
    written by permutations.py).
Fig. 11 (fig_validation_mechanisms): (a) shadowing shift of grating 2 behind grating 1 against the closed-form
    rule (eq:lawA), from the ruleA records of runs.py when present in records/, otherwise from
    cache/ruleA_shift.npz, (b) beat-limited SNR of the returns sharing one wavelength step, coherent model
    (markers, cache/beat_floor_rows.npy from beat_floor.py) against (eq:beat) (lines).
Output: ../figs/fig_validation_curvature and fig_validation_mechanisms (pdf, png).
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
import s70_field as F
import figstyle as FS

RECORDS = HERE / 'records'
CACHE = HERE / 'cache'
OUT = HERE.parent / 'figs'
BLUE, ORANGE = FS.BLUE, FS.VERM
FS.apply()
plt.rcParams.update({'pdf.fonttype': 42})


def save(fig, name):
    fig.savefig(OUT / (name + '.pdf'), bbox_inches='tight', pad_inches=.035)
    fig.savefig(OUT / (name + '.png'), bbox_inches='tight', pad_inches=.035)
    plt.close(fig)
    print('saved', name)


# ---------------------------------------------------------------- Fig. 11
fig, axs = plt.subplots(1, 2, figsize=(3.5, 1.6), layout='constrained')
# (a) shadowing shift: exact coefficient of the paper, no fitted gain
ax = axs[0]
u = np.linspace(-3.8, 3.8, 300)
ax.plot(u, -(4/3)*np.sqrt(2/3)*u*np.exp(-u*u/3), color='#333333', lw=1.4)
r = 0.10
if (RECORDS / 'ruleA_D3000_lw300.npz').exists():
    import runs as V
    z0 = np.load(RECORDS / 'ruleA_D3000_lw300.npz', allow_pickle=True)
    sig, c0 = V.width_of(z0, 1), V.centre_of(z0, 1)
    ds = np.array([-260, -130, 0, 50, 100, 130, 170, 200, 260, 330, 400])
    err = np.array([V.centre_of(np.load(RECORDS / ('ruleA_D%+04d_lw300.npz' % d), allow_pickle=True), 1) - c0 for d in ds])
else:
    zc = np.load(CACHE / 'ruleA_shift.npz')
    ds, err, sig = zc['detuning_pm'], zc['shift_pm'], float(zc['sigma_pm'])
    print('ruleA records absent, panel (a) from cache/ruleA_shift.npz')
ax.plot(ds/sig, err/(r*sig), 'o', color=ORANGE, ms=4, mfc=ORANGE)
pred = -(4/3)*np.sqrt(2/3)*r*ds*np.exp(-ds*ds/(3*sig*sig))
print('Rule A: sigma %.1f pm, max |shift| %.1f pm, RMSE against the rule %.2f pm' % (sig, abs(err).max(), np.sqrt(np.mean((err-pred)**2))))
ax.axhline(0, color='.6', lw=.6)
ax.set(xlabel=r'Difference $\Delta\lambda_{jk}/\sigma$', ylabel=r'Error $\delta\lambda_k/(R\sigma)$', xlim=(-3,4))
ax.legend(handles=[Line2D([],[],color='#333333',label='shadowing shift'),
    Line2D([],[],color=ORANGE,marker='o',ls='',label='$R_0=10\\%$')],
    loc='upper right', fontsize=5, handlelength=1.2, borderaxespad=0.3)
# (b) beat floor of the returns sharing one wavelength step: coherent model (markers) against (eq:beat) (lines)
rows = np.load(CACHE / 'beat_floor_rows.npy')          # columns: K, linewidth [Hz], chip rate [Hz], model floor, closed-form floor
N_JLT, R_JLT, EPS_JLT = 511, 0.01, 0.30
def beat_formula(K, dnu, B):
    P = (1 - R_JLT) ** (2 * np.arange(K))
    eta = (2 / np.pi) * np.arctan(B / (2 * dnu))
    return (1 + EPS_JLT) / (1 - EPS_JLT) * np.sqrt(eta * (P.sum() ** 2 - (P ** 2).sum()) / N_JLT)
b = axs[1]
dnu = np.logspace(np.log10(0.15e9), np.log10(5e9), 200)
for (K, B, mkr, lab), col in zip([(50, 31.25e6, 'o', '50 gratings, 31 Mchip/s'), (50, 100e6, 's', '50 gratings, 100 Mchip/s'), (17, 31.25e6, '^', '17 gratings, 31 Mchip/s')], [BLUE, ORANGE, FS.GREEN]):
    b.plot(dnu / 1e9, [20 * np.log10(1.0 / beat_formula(K, d, B)) for d in dnu], '-', color=col, lw=1.0)
    sel = (rows[:, 0] == K) & np.isclose(rows[:, 2], B)
    b.plot(rows[sel, 1] / 1e9, 20 * np.log10(1.0 / rows[sel, 3]), mkr, color=col, ms=3.6, mfc='white', mew=0.9, label=lab)
b.set_xscale('log')
b.set_xlim(0.15, 5); b.set_ylim(-2, 34)
b.set_xticks([0.3, 1, 3]); b.set_xticklabels(['0.3', '1', '3'])
b.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
b.set(xlabel='Source linewidth [GHz]', ylabel='Beat-limited SNR [dB]')
b.legend(fontsize=5, loc='upper left', handlelength=1.4, borderaxespad=0.3)
FS.letter(axs[0], 'a'); FS.letter(axs[1], 'b')
save(fig, 'fig_validation_mechanisms')

# ---------------------------------------------------------------- Fig. 3
fig = plt.figure(figsize=(3.5, 3.7), layout='constrained')
gs = fig.add_gridspec(2, 2, height_ratios=[1.9, 1.8])
ax = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
axc = fig.add_subplot(gs[1, :])
xx = np.linspace(-350,350,1401)
# (a) direct return of the 8th co-tuned grating, coupled-mode spectra (solid) and Gaussian lines (dotted)
for r, col, lsty in [(.04, BLUE, (0, (5, 2))), (1/15, FS.ORANGE, (0, (4, 1.5, 1, 1.5))), (.1, ORANGE, '-')]:
    length, kappa, _, _ = F.calibrate('bl250', r)
    q = abs(F.grating_r(F.LAM0+xx*1e-12, F.LAM0, length, kappa, sections=160))**2
    direct = q*(1-q)**14
    gauss=r*np.exp(-.5*(xx/(250/2.354820045))**2)
    gd=gauss*(1-gauss)**14
    ax[0].plot(xx, direct/direct.max(), color=col, lw=1.4, ls=lsty, label=f'$R_0={r*100:.2g}\\%$')
    ax[0].plot(xx, gd/gd.max(), ls=(0, (1.2, 1.6)), color=col, lw=.9)
FS.letter(ax[0], 'a')
ax[0].set(xlabel='Offset [pm]', ylabel='Normalized $A_8$', xlim=(-300,300), ylim=(-.03,1.15))
ax[0].legend(loc='lower center', ncol=1, fontsize=6, handlelength=2.6, borderaxespad=0.3)
ax[0].text(.03,.97, r'$R_c=1/15$', transform=ax[0].transAxes, va='top', fontsize=6)
# (b) the same eight gratings of unequal width, narrow first and wide first
for tag, col, widths in [('narrowfirst',ORANGE,np.linspace(175,325,8)), ('widefirst',BLUE,np.linspace(325,175,8))]:
    length,kappa,_,_=F.calibrate('bl250',.1)
    q=[]
    for w in widths:
        q.append(abs(F.grating_r(F.LAM0+xx*1e-12,F.LAM0,length*250/w,kappa*w/250,sections=160))**2)
    direct=q[-1]*np.prod((1-np.asarray(q[:-1]))**2,axis=0)
    lab='Narrow first' if tag=='narrowfirst' else 'Wide first'
    qg=.1*np.exp(-.5*(xx[None,:]/(widths[:,None]/2.354820045))**2)
    gd=qg[-1]*np.prod((1-qg[:-1])**2,axis=0)
    ic = np.argmin(np.abs(xx))          # common value at the Bragg wavelength: the product of transmissions there does not depend on the order
    ax[1].plot(xx,direct/direct[ic],color=col,lw=1.4,ls=('-' if tag=='narrowfirst' else (0, (5, 2))),label=lab)
    ax[1].plot(xx,gd/gd[ic],ls=(0, (1.2, 1.6)),color=col,lw=.9)
    gamma=1/widths**2
    omega=np.r_[0,np.cumsum(2*.1/.9*gamma)[:-1]]/gamma
    print('%s: Omega_max %.2f' % (tag, omega.max()))
FS.letter(ax[1], 'b')
ax[1].set(xlabel='Offset [pm]', ylabel='$A_8/A_8(\\lambda_B)$', xlim=(-300,300), ylim=(-.03,1.85))
ax[1].text(150, 1.72, r'$\Omega_{\max}=3.17$', color=ORANGE, fontsize=6, ha='center', va='center')
ax[1].text(0, 1.12, r'$\Omega_{\max}=0.76$', color=BLUE, fontsize=6, ha='center', va='bottom')
ax[1].legend(loc='lower center', fontsize=6, handlelength=2.6, borderaxespad=0.3)
ax[1].text(.03,.97,'$R_0=10\\%$',transform=ax[1].transAxes,va='top',fontsize=6)
# (c) all 40320 orders of the same eight gratings (permutations.py)
zp = np.load(CACHE / 'order_perm.npz')
perms, om_all, nsplit, off_all = zp['perms'], zp['omega_max'], zp['n_split'], zp['shift_max']
i_wide = int(np.where((perms == np.arange(8)[::-1]).all(axis=1))[0][0])
i_narrow = int(np.where((perms == np.arange(8)).all(axis=1))[0][0])
single = nsplit == 0
sel = np.random.default_rng(1).permutation(len(perms))
axc.scatter(om_all[sel][~single[sel]], off_all[sel][~single[sel]], s=3, c=ORANGE, lw=0, alpha=.3, rasterized=True)
axc.scatter(om_all[sel][single[sel]], off_all[sel][single[sel]], s=4, c=BLUE, lw=0, alpha=.9, rasterized=True, zorder=4)
omg = np.linspace(0.3, 0.97, 200)
axc.plot(omg, 5 * omg / (1 - omg), 'k--', lw=.9)
axc.axvline(1.0, color='0.5', ls=':', lw=.8)
axc.plot(om_all[i_wide], off_all[i_wide], 'o', ms=4.5, mfc='white', mec='k', mew=.8, zorder=6)
axc.plot(om_all[i_narrow], off_all[i_narrow], 's', ms=4.5, mfc='white', mec='k', mew=.8, zorder=6)
axc.set_yscale('log')
axc.set(xlabel=r'$\Omega_{\max}$ of the order', ylabel='Offset of maximum [pm]', xlim=(0.6, 3.4), ylim=(2, 3000))
axc.legend(handles=[Line2D([], [], ls='', marker='o', color=ORANGE, ms=3.2, label='Split return (%d orders)' % (~single).sum()),
                    Line2D([], [], ls='', marker='o', color=BLUE, ms=3.2, label='Single peak (%d orders)' % single.sum()),
                    Line2D([], [], color='k', ls='--', lw=.9, label=r'$5\,\Omega_{\max}/(1-\Omega_{\max})$'),
                    Line2D([], [], ls='', marker='o', ms=4.5, mfc='white', mec='k', mew=.8, label='wide first, as in (b)'),
                    Line2D([], [], ls='', marker='s', ms=4.5, mfc='white', mec='k', mew=.8, label='narrow first, as in (b)')],
           fontsize=5.5, loc='lower right', handlelength=1.8, borderaxespad=.3)
# insets: return of the grating with the largest Omega_k for three orders, normalized to the value at lambda_B
w8 = np.linspace(175, 325, 8)
length10, kappa10, _, _ = F.calibrate('bl250', .1)
qw = {w: abs(F.grating_r(F.LAM0 + xx * 1e-12, F.LAM0, length10 * 250 / w, kappa10 * w / 250, sections=160)) ** 2 for w in w8}
for target, pos in ((0.95, (0.27, 0.64)), (1.5, (0.48, 0.64)), (2.5, (0.69, 0.64))):
    i_sel = int(np.argmin(np.abs(om_all - target)))
    order = w8[perms[i_sel]]
    omk = np.r_[0, np.cumsum(2 * .1 / .9 * (1 / order ** 2))[:-1]] / (1 / order ** 2)
    k = int(np.argmax(omk))
    ret = qw[order[k]] * np.prod([(1 - qw[w]) ** 2 for w in order[:k]], axis=0)
    ret = ret / ret[np.argmin(np.abs(xx))]
    ins = axc.inset_axes([pos[0], pos[1], 0.19, 0.34])
    ins.plot(xx, ret, color=ORANGE if om_all[i_sel] >= 1 else BLUE, lw=1.0)
    ins.set(xlim=(-300, 300), ylim=(0, 1.75), xticks=[], yticks=[])
    ins.text(.5, .95, r'$\Omega_{\max}=%.2g$' % om_all[i_sel], transform=ins.transAxes, ha='center', va='top', fontsize=5)
    axc.plot(om_all[i_sel], off_all[i_sel], 'D', ms=3.5, mfc='white', mec='k', mew=.7, zorder=7)
    axc.annotate('', xy=(om_all[i_sel], off_all[i_sel]), xytext=(pos[0] + 0.095, pos[1]), textcoords='axes fraction',
                 arrowprops=dict(arrowstyle='-', color='0.3', lw=.6, shrinkA=0, shrinkB=2.5))
FS.letter(axc, 'c')
print('orders: %d, single peak %d, Omega_max < 1: %d, wide first Omega %.2f offset %.1f pm, narrow first Omega %.2f offset %.1f pm' % (
    len(perms), single.sum(), (om_all < 1).sum(), om_all[i_wide], off_all[i_wide], om_all[i_narrow], off_all[i_narrow]))
save(fig, 'fig_validation_curvature')
