# python_selfcal_cdm — samokalibrująca interrogacja CDM (badanie symulacyjne)

Lekkie, reprodukowalne symulacje w Pythonie do artykułu o **tłumieniu w trybie
wspólnym błędu FM→AM od chirpu źródła** w interrogacji CDM sieci siatek Bragga z
użyciem **współ-kodowanych siatek referencyjnych**. Uzupełnienie pełnego
symulatora MATLAB w katalogu nadrzędnym (`../`).

## Wymagania
- Python 3.10+, `numpy`, `scipy`, `matplotlib`
- `pip install numpy scipy matplotlib`

## Pliki
| Plik | Opis |
|------|------|
| `common.py` | wspólna fizyka: kody Gold/Kasami (konstrukcja decymacyjna), widma FBG (Gauss / tanh²), model chirpu i konwersji FM→AM, estymatory piku (Gauss / centroid), nieliniowy przeciąg VCSEL, linijka MZI (k-clock) |
| `s1_codes.py` | właściwości korelacyjne kodów; podłoga listków bocznych vs długość kodu |
| `s2_despread.py` | rozplot CDM nakładających się widmowo siatek + kanały referencyjne |
| `s3_residual.py` | **wynik główny**: rezyduum błędu λ_B po korekcie referencjami vs chirp/FWHM i niedopasowanie nachylenia |
| `s4_baselines.py` | porównanie strategii kalibracji: brak / MZI k-clock / referencje / referencje+MZI |
| `s5_ghosts.py` | odporność kotwicy referencyjnej na wielokrotne odbicia (ghosty) |

## Uruchomienie
```bash
python3 s1_codes.py     # i analogicznie s2..s5; figury zapisywane do figs/
```

## Kluczowe wyniki (parametry domyślne)
- Kody Gold N=127: maks. korelacja wzajemna **0,134 = granica teoretyczna 17/127**; Kasami N=63: 9/63.
- Rozplot CDM (idealny tor, bez chirpu): MAE ≈ **9 pm**.
- Korekta referencjami: surowy błąd ≈ **50 pm** @Δ=0,5·FWHM → **3,7 pm** (4 referencje). Rezyduum rośnie z chirp/FWHM i maleje z liczbą referencji; podłoga = niedopasowanie nachylenia widm.
- Strategie (średnia z 40 realizacji): brak **110 pm**, tylko MZI k-clock **50 pm** (zostaje FM→AM), tylko referencje **39 pm**, referencje+MZI **3,8 pm**.
- Ghosty: pozycja kotwicy (dopasowanie Gaussa) praktycznie niezmiennicza (**~0 pm @R=9%**), naiwny centroid pociągany (**35 pm**).

Jednostki: częstotliwość optyczna w GHz, przeliczenie 1 GHz ≈ 8 pm przy 1550 nm.
