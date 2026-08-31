"""s45_lineshape.py - the Gaussian line model of the paper against the
reflectivity of a uniform grating from coupled-mode theory (Erdogan 1997).

Both lines share the peak R0 = 10 % and the FWHM = 250 pm. For the uniform
grating that fixes kappa L = atanh(sqrt(R0)) = 0.33 and, at 1545 nm with
n_eff = 1.45, a length of 3.0 mm. The symbols of the line model are marked on
the Gaussian: R0 at the peak, lambda_B at its center, the FWHM at half height
and sigma at R0 exp(-1/2). The side lobes of the uniform grating are what
apodization removes in gratings made for sensing, so the Gaussian stands for
the apodized line.

Output: figs/fig_s45_lineshape.pdf and .png, drawn at the printed width
(0.85 of a 3.5-in IEEE column) in the figstyle look.
"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import figstyle as FS
FS.apply()

R0 = 0.10
Q = np.arctanh(np.sqrt(R0))            # kappa L
FWHM = 250.0                           # pm
SIG = FWHM / (2 * np.sqrt(2 * np.log(2)))
LAM = 1545e-9
NEFF = 1.45


def erdogan(x, q):
    """Uniform-grating reflectivity, x = sigma_hat * L, q = kappa L."""
    x = np.asarray(x, float)
    r = np.empty_like(x)
    a = np.abs(x) < q
    s = np.sqrt(q**2 - x[a]**2)
    r[a] = np.sinh(s)**2 / (np.cosh(s)**2 - x[a]**2 / q**2)
    b = ~a
    s = np.sqrt(x[b]**2 - q**2)
    r[b] = np.sin(s)**2 / (x[b]**2 / q**2 - np.cos(s)**2)
    return r


# half-maximum detuning in x units, bisection between the peak and the first zero
lo, hi = Q, np.sqrt(Q**2 + np.pi**2)
for _ in range(80):
    mid = 0.5 * (lo + hi)
    if erdogan([mid], Q)[0] > 0.5 * R0:
        lo = mid
    else:
        hi = mid
XH = 0.5 * (lo + hi)
L = XH * LAM**2 / (2 * np.pi * NEFF * (FWHM / 2) * 1e-12)
print('kappa L = %.3f, L = %.2f mm, kappa = %.0f 1/m' % (Q, L * 1e3, Q / L))

dl = np.linspace(-600, 600, 2401)
r_uni = 100 * erdogan(dl * XH / (FWHM / 2), Q)
r_gau = 100 * R0 * np.exp(-0.5 * (dl / SIG)**2)
side = r_uni[dl > 300].max() / (100 * R0) * 100
print('side lobe %.1f %% of the peak' % side)

c_uni, c_gau, c_mark = FS.VERM, FS.BLUE, FS.GREY
fig, ax = plt.subplots(figsize=(2.98, 1.95))

# markers for R0 and lambda_B
ax.plot([-600, 0], [10, 10], color=c_mark, lw=0.6, ls=(0, (3, 2)), zorder=1)
ax.plot([0, 0], [0, 10], color=c_mark, lw=0.6, ls=(0, (3, 2)), zorder=1)
ax.plot(dl, r_uni, color=c_uni, label='uniform grating')
ax.plot(dl, r_gau, color=c_gau, label='Gaussian model')

ax.text(-585, 10.2, r'$R_{0,k}$', color=c_mark, ha='left', va='bottom', fontsize=7.5)
ax.text(-14, 2.3, r'$\lambda_{B,k}$', color=c_mark, ha='right', va='center', fontsize=7.5)

# FWHM at half height, sigma at R0 exp(-1/2)
half = FWHM / 2
ax.plot([-half, half], [5, 5], color=c_gau, lw=0.8, zorder=4)
for x in (-half, half):
    ax.plot([x, x], [4.5, 5.5], color=c_gau, lw=0.8, zorder=4)
ax.text(half + 18, 5, r'$\mathrm{FWHM}_k$', color=c_gau, ha='left', va='center',
        fontsize=7.5, bbox=dict(fc='white', ec='none', alpha=0.85, pad=0.6))
ys = 10 * np.exp(-0.5)
ax.plot([0, SIG], [ys, ys], color=c_gau, lw=0.8, zorder=4)
ax.plot([SIG, SIG], [ys - 0.5, ys + 0.5], color=c_gau, lw=0.8, zorder=4)
ax.text(SIG + 18, ys + 0.55, r'$\sigma_k$', color=c_gau, ha='left', va='center',
        fontsize=7.5, bbox=dict(fc='white', ec='none', alpha=0.85, pad=0.6))

ax.set_xlim(-600, 600)
ax.set_ylim(0, 13.8)
ax.set_xticks([-400, -200, 0, 200, 400])
ax.set_yticks([0, 5, 10])
ax.set_xlabel(r'$\lambda-\lambda_{B,k}$ (pm)')
ax.set_ylabel(r'$R_k(\lambda)$ (%)')
ax.legend(loc='upper right', handlelength=1.6, borderaxespad=0.4,
          labelspacing=0.25, handletextpad=0.5)

fig.tight_layout(pad=0.15)
fig.savefig('figs/fig_s45_lineshape.png', dpi=300, bbox_inches='tight', pad_inches=0)
fig.savefig('figs/fig_s45_lineshape.pdf', bbox_inches='tight', pad_inches=0)
print('saved figs/fig_s45_lineshape.pdf')
