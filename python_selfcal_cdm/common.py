"""
common.py - wspolna fizyka do badania symulacyjnego: samokalibrujaca sie
interrogacja CDM sieci siatek Bragga z bezposrednio modulowanym, przestrajanym
VCSEL-em. Projekt PN/01/0321/2022 (ZPSSC).

Jednostki: czestotliwosc optyczna w GHz wzgledem srodka pasma; przesuniecia
dlugosci fali w pm, przelicznik 1 GHz ~ 8 pm przy 1550 nm.
"""
import numpy as np
from scipy.signal import max_len_seq
from scipy.optimize import curve_fit
import warnings; warnings.filterwarnings("ignore")

GHZ_PER_PM = 0.125            # 1 pm -> 0.125 GHz @1550 nm  (dnu = c*dlam/lam^2)
PM_PER_GHZ = 1.0 / GHZ_PER_PM # ~8.0 pm na 1 GHz

# ----------------------------------------------------------------------------
# Sekwencje kodowe (m-ciagi, Gold, Kasami) - konstrukcja decymacyjna
# ----------------------------------------------------------------------------
def _mls01(nbits):
    seq, _ = max_len_seq(nbits)      # 0/1, dlugosc N=2^n-1
    return seq.astype(int)

def _decimate01(seq, q):
    N = len(seq)
    idx = (q * np.arange(N)) % N
    return seq[idx]

def _to_pm1(seq01):
    return 1.0 - 2.0 * seq01         # 0->+1, 1->-1  (bipolarny)

def gold_set(nbits):
    """Rodzina Gold dla nieparzystego n (para preferowana przez decymacje q=2^((n+1)/2)+1)."""
    N = 2**nbits - 1
    u = _mls01(nbits)
    q = 2**((nbits + 1)//2) + 1
    v = _decimate01(u, q)
    codes = [u.copy(), v.copy()]
    for k in range(N):
        codes.append((u ^ np.roll(v, k)))
    return np.array([_to_pm1(c) for c in codes])

def kasami_small(nbits):
    """Maly zbior Kasamiego dla parzystego n (decymacja q=2^(n/2)+1)."""
    assert nbits % 2 == 0
    N = 2**nbits - 1
    u = _mls01(nbits)
    q = 2**(nbits//2) + 1
    w = _decimate01(u, q)            # okres 2^(n/2)-1, powielony do N
    codes = [u.copy()]
    for k in range(2**(nbits//2) - 1):
        codes.append(u ^ np.roll(w, k))
    return np.array([_to_pm1(c) for c in codes])

def periodic_xcorr(a, b):
    """Cykliczna korelacja wzajemna (znormalizowana do dlugosci)."""
    n = len(a)
    A = np.fft.fft(a); B = np.fft.fft(b)
    return np.fft.ifft(A * np.conj(B)).real / n

def code_corr_stats(codes):
    """Szczytowa autokorelacja boczna i szczytowa korelacja wzajemna (znorm. do N)."""
    M, N = codes.shape
    auto_side = 0.0; cross = 0.0
    for i in range(M):
        c = periodic_xcorr(codes[i], codes[i])
        auto_side = max(auto_side, np.max(np.abs(c[1:])))   # poza pikiem w 0
        for j in range(i+1, M):
            cc = periodic_xcorr(codes[i], codes[j])
            cross = max(cross, np.max(np.abs(cc)))
    return auto_side, cross   # piki w 0 sa rowne 1 (=N/N)

# ----------------------------------------------------------------------------
# Widma siatek Bragga
# ----------------------------------------------------------------------------
def fbg_gauss(nu, nu_b, fwhm):
    sig = fwhm / 2.35482
    return np.exp(-0.5 * ((nu - nu_b) / sig)**2)

def fbg_tanh(nu, nu_b, fwhm, n_side=0.0):
    """Przyblizenie tanh^2 (coupled-mode) z opcjonalna asymetria n_side.
    Asymetria modeluje rzeczywista, lekko niesymetryczna charakterystyke siatki."""
    sig = fwhm / 2.35482
    x = (nu - nu_b) / sig
    base = 1.0 / np.cosh(x)**2
    if n_side != 0.0:
        base = base * (1.0 + n_side * np.tanh(x))   # lekka asymetria zboczy
        base = np.clip(base, 0.0, None)
    return base / base.max()

# ----------------------------------------------------------------------------
# Chirp / konwersja FM->AM
# ----------------------------------------------------------------------------
def chirp_kernel(delta, skew=0.0, npts=401, span=4.0):
    """Rozklad chwilowych odchylek czestotliwosci podczas 'zapalonych' chipow kodu.
    delta - rozmiar (odchylenie std) ekskursji chirpu [GHz]; skew - skosnosc
    (modulacja bezposrednia daje niesymetryczny chirp transient vs adiabatyczny)."""
    if delta <= 0:
        d = np.array([0.0]); p = np.array([1.0]); return d, p
    d = np.linspace(-span*delta, span*delta, npts)
    p = np.exp(-0.5 * (d/delta)**2)
    if skew != 0.0:
        from scipy.special import erf
        p = p * (1.0 + erf(skew * d / (np.sqrt(2)*delta)))   # rozklad skosny
    p = np.clip(p, 0, None); p /= p.sum()
    return d, p

def fmam_readout(nu_grid, nu_b, fwhm, delta, mean_off=0.0, skew=0.0,
                 shape='tanh', asym=0.0):
    """Efektywny, zdespreadowany odczyt widmowy siatki z rozmyciem FM->AM.
    Chirp powoduje, ze przy pozycji przeciagu nu_s odbiornik usrednia widmo
    siatki po rozkladzie chwilowych czestotliwosci (kernel chirpu)."""
    d, p = chirp_kernel(delta, skew=skew)
    d = d + mean_off
    if shape == 'gauss':
        spec = lambda x: fbg_gauss(x, nu_b, fwhm)
    else:
        spec = lambda x: fbg_tanh(x, nu_b, fwhm, n_side=asym)
    out = np.zeros_like(nu_grid)
    for di, pi in zip(d, p):
        out += pi * spec(nu_grid + di)
    return out

# ----------------------------------------------------------------------------
# Estymatory polozenia piku
# ----------------------------------------------------------------------------
def _gauss(x, a, mu, sig, c):
    return a * np.exp(-0.5*((x-mu)/sig)**2) + c

def gauss_fit_peak(x, y, thr=0.3, win_frac=0.14):
    """Dopasowanie Gaussa w lokalnym oknie wokol piku (odporne na listki/cross)."""
    x=np.asarray(x,float); y=np.asarray(y,float)
    i0=int(np.argmax(y)); W=max(4,int(len(x)*win_frac))
    lo=max(0,i0-W); hi=min(len(x),i0+W+1)
    xs=x[lo:hi]; ys=y[lo:hi]
    ym=ys.max(); base=float(np.min(ys))
    m=ys>=base+thr*(ym-base)
    if m.sum()<4: return centroid(xs,ys,thr)
    p0=[ym-base, xs[int(np.argmax(ys))], (xs[-1]-xs[0])/6.0, base]
    try:
        lob=[0.0,xs[0],1e-6,-np.inf]; hib=[np.inf,xs[-1],xs[-1]-xs[0],np.inf]
        popt,_=curve_fit(_gauss,xs[m],ys[m],p0=p0,bounds=(lob,hib),maxfev=20000)
        return float(popt[1])
    except Exception:
        return centroid(xs,ys,thr)

def centroid(x, y, thr=0.3):
    ym = y.max(); m = y >= thr*ym
    w = y[m] - thr*ym
    return np.sum(x[m]*w)/np.sum(w)

# ----------------------------------------------------------------------------
# Nieliniowy przeciag VCSEL i interferometr MZI (k-clock)
# ----------------------------------------------------------------------------
def vcsel_sweep(v, span_ghz, quad=0.30, cube=0.0):
    """nu(V) ~ kwadratowo-szescienna funkcja napiecia strojenia HCG-VCSEL
    (BW10: Delta lambda ~ kwadratowo z napieciem strojenia). v w [0,1]."""
    base = v + quad*(v**2 - v) + cube*(v**3 - v)
    base = (base - base.min())/(base.max()-base.min())
    return base * span_ghz

def mzi_kclock_estimate(nu_true, fsr_ghz):
    """Linijka czestotliwosci MZI: zwraca 'znaczniki' rownych przyrostow nu
    (zera prazkow). Idealna linijka czestotliwosci zrodla (lacznie z chirpem
    wspolnym), ale BEZ informacji o FM->AM na zboczu siatki (efekt detekcji)."""
    n0 = nu_true.min(); n1 = nu_true.max()
    ticks = np.arange(n0, n1, fsr_ghz)
    return ticks
