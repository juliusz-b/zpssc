"""Schematy pogladowe: (Fig.1) diagram blokowy systemu, (Fig.2) odwzorowanie
dlugosc fali - czas dla kodowanego przeciagu."""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
# ---- Fig.1 block diagram ----
fig,ax=plt.subplots(figsize=(11,3.2)); ax.axis('off'); ax.set_xlim(0,11); ax.set_ylim(0,3)
def box(x,y,w,h,t,fc='#eaf1fb'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.04,rounding_size=0.08",fc=fc,ec='#1f4e79',lw=1.4))
    ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=8.5)
def arrow(x1,y1,x2,y2):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=12,lw=1.3,color='#333'))
box(0.1,1.1,1.8,0.9,"Przestrajalny\nVCSEL 1550 nm\n+ kod (mod. natez.)")
box(2.3,1.2,1.0,0.7,"Cyrkulator")
box(3.7,0.5,3.2,2.0,"Lancuch FBG\n S1  R1  S2  R2  S3\n(czujnikowe + referencyjne\nna modulach Peltiera)",fc='#fdf3e7')
box(7.3,1.2,1.0,0.7,"APD\nInGaAs")
box(8.7,0.9,2.1,1.3,"STM32H7:\nrozplot (FFT)\n+ kalibracja ref.\n+ ekstrakcja lambda",fc='#eafaf0')
arrow(1.9,1.55,2.3,1.55); arrow(3.3,1.55,3.7,1.55)
arrow(5.3,0.5,5.3,0.2); arrow(5.3,0.2,7.8,0.2); arrow(7.8,0.2,7.8,1.2)  # reflections to APD
arrow(8.3,1.55,8.7,1.55)
ax.text(5.3,2.62,"Fig. 1. Schemat blokowy: kodowany przestrajalny VCSEL -> lancuch FBG (czujniki + referencje) -> APD -> rozplot i samokalibracja na wezle brzegowym.",ha='center',fontsize=7.5)
plt.tight_layout(); plt.savefig('figs/fig0_blockdiagram.png',dpi=150,bbox_inches='tight')
# ---- Fig.2 wavelength-time mapping ----
fig,ax=plt.subplots(figsize=(8,3.2))
t=np.linspace(0,1,1000); nu=2.0*t-1.0  # ramp sweep (normalized band)
ax.plot(t,nu,'b',lw=1.6,label='przeciag VCSEL $\\nu_s(t)$ (nieliniowy)')
nu_nl=2.0*t-1.0+0.18*np.sin(2*np.pi*t)  # nonlinear
ax.plot(t,nu_nl,'b--',lw=1.0,alpha=0.7,label='rzeczywisty (nieliniowy)')
for lb,name in [(-0.5,'S1'),(0.0,'R1'),(0.45,'S2')]:
    tc=(lb+1)/2.0
    ax.axhline(lb,color='0.85',lw=0.6)
    ax.plot([tc],[lb],'r^',ms=8); ax.text(tc,lb+0.06,name+' ($\\lambda_B$)',fontsize=8)
# code modulation cartoon on top
tt=np.linspace(0.02,0.18,200); code=0.9+0.06*((np.sin(2*np.pi*tt*60)>0)*2-1)
ax.plot(tt,code,'k',lw=0.7); ax.text(0.1,0.99,'kod c(t)',fontsize=7,ha='center')
ax.set_xlabel('czas przeciagu (znorm.)'); ax.set_ylabel('odstrojenie optyczne (znorm.)')
ax.set_title('Fig. 2. Odwzorowanie dlugosc fali - czas: kazda siatka odbija, gdy $\\nu_s(t)$ trafia w jej rezonans;\nkod nadaje zysk rozplotu, a nieliniowosc przeciagu i chirp wymagaja kalibracji.',fontsize=7.5)
ax.legend(fontsize=7,loc='lower right')
plt.tight_layout(); plt.savefig('figs/fig0_wavelength_time.png',dpi=150,bbox_inches='tight')
print("saved fig0_blockdiagram.png, fig0_wavelength_time.png")
