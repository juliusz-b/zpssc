# python_selfcal_cdm — samokalibrująca interrogacja CDM (badanie symulacyjne)

Lekkie, reprodukowalne symulacje w Pythonie do artykułu JCR (EN) o samokalibrującej
się interrogacji CDM sieci siatek Bragga z tanim, bezpośrednio modulowanym,
przestrajanym VCSEL-em. Uzupełnienie symulatora MATLAB (`../`). Kod i opisy po angielsku.

> **STAN:** po rundzie 6 recenzji (major revision), w rewizji. Handoff i decyzje:
> `../../Realizacja zadań/Artykul_JCR_2026/STAN_ARTYKULU_2026-06.md` oraz `RECENZJE_I_DECYZJE_2026.md`.

## Wymagania
- Python 3.10+, `numpy`, `scipy`, `matplotlib` (`pip install numpy scipy matplotlib`)

## Pliki
| Plik | Opis |
|------|------|
| `common.py` | wspólna fizyka + parametry sprzętu: kody Gold/Kasami, widma FBG (Gauss/tanh²), model chirpu i FM→AM, estymatory piku (Gauss/centroid), nieliniowy przeciąg VCSEL, linijka MZI |
| `s0_schematics.py` | schematy poglądowe: diagram blokowy + odwzorowanie długość fali-czas |
| `s1_codes.py` | właściwości korelacyjne kodów; podłoga listków vs długość kodu (Fig.3) |
| `s2_despread.py` | rozplot CDM nakładających się siatek + kanały referencyjne (Fig.4, waterfall) |
| `s3_residual.py` | rezyduum błędu λ_B po korekcie referencjami vs chirp/FWHM i niedopasowanie nachylenia (Fig.5) |
| `s4_baselines.py` | strategie kalibracji: brak / MZI / referencje / referencje+MZI (Fig.6) |
| `s5_ghosts.py` | odporność kotwicy na wielokrotne odbicia (Fig.7) |
| `s6_experiment.py` | predykcja makiety (FBGS DTG-A3A4, 1545 nm, 250 pm, R=10%) (Fig.8) |
| `s7_axis_nonlinearity.py` | kalibracja nieliniowej osi λ(V) HCG; szerokie pasmo vs wąska makieta (Fig.9) |
| `s_setup_figure.py` | schemat stanowiska eksperymentalnego (Fig.10) |
| `s8_noise.py` | podłoga szumowa: jitter centroidu vs SNR i reflektancja |
| `s9_coherent.py` | długość koherencji z linii CW; przypadek koherentny vs niekoherentny |
| `s10_control.py` | **kontrola nowości**: referencja ko-kodowana vs zwykła (decydujące) |
| `s11_sensitivity.py` | wrażliwość rezyduum na rozrzut niesymetrii i kształt (tanh² vs sinc²) |

## Uruchomienie
```bash
python3 s1_codes.py     # analogicznie s2..s11; figury -> figs/
```

## Kluczowe wyniki
- Kody Gold N=127: korelacja wzajemna **0,134 = granica 17/127**; podłoga maleje z długością kodu.
- Rozplot CDM (idealny tor): MAE ≈ **9 pm**.
- Referencje usuwają wspólny offset chirpu → **podłoga rezyduum ~kilka pm** (niedopasowanie nachylenia). Surowy błąd FM-AM **skaluje się z chirpem** (chirp do zmierzenia na OSA; NIE zakładać konkretnej wartości).
- Oś λ(V) (s7): szerokie pasmo wymaga 3 ref lub MZI; wąska makieta — 2 ref wystarczą (~9 pm).
- Ghosty: kotwica (Gauss) niezmiennicza; koherencja (s9): L_coh ≈ 0,22 m ≪ 4 m → **niekoherentnie**.
- Szum (s8): przy SNR ≥16 dB systematyka rządzi.
- **Kontrola nowości (s10): ko-kodowanie referencji to zaleta SYSTEMOWA (bez dodatkowego sprzętu), NIE poprawa dokładności** — przy innych kanałach nawet szkodzi (cross-talk). Zmienia tezę pracy.

Jednostki: częstotliwość optyczna w GHz, 1 GHz ≈ 8 pm przy 1550 nm.
