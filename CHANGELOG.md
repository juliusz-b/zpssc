# Changelog - zpssc simulator

Format: data | zadanie | pliki | opis zmian

---

## 2026-01-15 | Optymalizacja Gaussian fit - 21.6 pm MAE

### Wyniki
- N_s=5: Default 391 → Opt Centroid 56 → **Opt Gaussian 21.6 pm (−94.5%)**
- N_s=10: Default 402 → Opt Centroid 93 → **Opt Gaussian 49.7 pm (−87.7%)**
- Gaussian fit preferuje wyższy spread (0.61nm vs 0.34nm) i niższe R (2.6%)
- **21.6 pm - najlepszy wynik w projekcie, sub-25pm dokładność**

### Pliki
- `results/WP3/optimization/gauss_optimization_results.mat`
- `results/WP3/optimization/gauss_vs_centroid_optimization.fig/.png`

---

## 2026-01-14 | Gaussian fit integracja + TDM per-siatka analiza

### Nowe pliki
- `src/system/runSimulationGauss.m` - runSimulation z Gaussian fit (spline 10× + lsqcurvefit)
- `src/system/optimObjectiveGauss.m` - funkcja celu dla GA/PSO z Gaussian fit
- `scripts/Zadanie4_GaussOptimization.m` - optymalizacja PSO z Gaussian fit (do uruchomienia ~5min)

### Wyniki
- Gaussian fit: p=8 −63%, p=10 −68% poprawy vs centroid
- Najlepszy wynik: p=10 + Gauss = **45.1 pm MAE**
- TDM per-siatka: ghosty zmieniają błędy ±400 pm per siatka (vs ±2 pm w CDM)
- CDM stabilizuje detekcję per siatka - nowy argument za CDM

---

## 2026-01-12 | Dokumentacja końcowa

### Podsumowanie prac grudniowych/styczniowych (~6.5h pracy)
- **14 skryptów MATLAB** w scripts/
- **28 plików wyników** (png, fig, mat) w results/
- **3 nowe funkcje core** (runSimulation, defaultParams, optimObjective)
- **Kluczowe odkrycia**: floor=sidelobes, CDM tłumi ghosty, Gaussian fit −88%, spectral spread jako parametr

---

## 2026-01-10 | FFE: finalny wniosek - CDM tłumi ghosty naturalnie

### Odkrycie
Ghosty R³≈17% (sygnał czasowy) są tłumione do **0.1–1.6%** po korelacji CDM. Kody Kasami są quasi-ortogonalne → korelacja krzyżowa ghost(kod_k) × kod_k' ≈ 0. CDM pełni rolę **naturalnego filtra anty-ghostingowego**.

FFE testowane: sample-rate (0%), sparse (0%), korelogram (6.9%), kanałowe kody1 (0%). Żadne nie daje istotnej poprawy bo ghosty już stłumione przez CDM. FFE miałoby sens w TDM (bez kodów).

Gaussian fit: **63–88% poprawy** MAE - jedyna skuteczna metoda zaawansowanego przetwarzania.

Wniosek dla projektu: **CDM >> TDM w kontekście anti-ghosting** - nowy argument za multipleksacją kodową.

---

## 2026-01-08 | Odtwarzanie widma FBG + testy FFE

### Nowe pliki
- `results/WP2/benchmark/spectrum_reconstruction_p8_vs_p10.fig/.png` - porównanie jakości widma
- `scripts/Zadanie_DFE_test.m` - testy SIC/DFE/FFE (zaktualizowany)

### Wyniki
- Jakość widma: p=8 peak/floor=3.5:1, p=10 peak/floor=5.1:1
- FFE RLS: 4 podejścia przetestowane, żadne nie daje stabilnej poprawy MAE - wymaga dalszych badań
- Gaussian fit nadal najskuteczniejszy: p=8→144pm (−63%), p=10→45pm (−88%)

---

## 2026-01-05 | Analiza SIC/DFE + Gaussian fit (PRZEŁOM)

### Wyniki
- SIC/DFE/LMS: nieskuteczne (0–2.4% poprawy) - błąd MAE to bias centroidu, nie crosstalk
- Więcej lambdas = gorzej (8λ: 145pm vs 64λ: 450pm) - więcej MAI
- **Gaussian fit + spline interpolation**: MAE p=8: 391→**144 pm (−63%)**, p=10: 139→**45 pm (−88%)**
- Najlepszy wynik w projekcie: **45.1 pm MAE** (p=10 + Gaussian fit)

### Nowe pliki
- `scripts/Zadanie_DFE_test.m` - testy SIC/DFE (negatywne, ale ważne wyniki)

---

## 2025-12-20 | Zadanie 5 (WP3.4) - porównanie przed/po optymalizacji

### Nowe pliki
- `scripts/Zadanie5_Comparison.m` - porównanie default vs GA/PSO optimum, 4 warianty
- `results/WP3/comparison/comparison_v1..v4.fig/.png` - szczegółowe wykresy (6 subplotów)
- `results/WP3/comparison/comparison_summary.fig/.png` - wykres zbiorczy
- `results/WP3/comparison/comparison_results.mat` - dane

### Wyniki
- N_s=5: 391→**56 pm** (−86%, PSO), R_opt=10.2%, spread=0.34nm
- N_s=10: 402→**93 pm** (−77%, GA), R_opt=0.4%, spread=0.14nm
- N_s=10+grad: 463→**130 pm** (−72%, PSO), R_opt=0.6%, spread=0 (gradient = naturalny spread)
- N_s=20: 396→**287 pm** (−28%, GA), R_opt=1.0% - p=10 nie wystarcza, potrzeba p=12

Zaskakujące: gradient "wyłącza" spectral spread; Nb=16 tylko bez gradientu; N_s=20 wymaga dłuższych kodów.

---

## 2025-12-18 | Zadanie 4 - wyniki pełnej optymalizacji GA vs PSO

### Wyniki (4 warianty, 5 zmiennych)
- N_s=5: PSO **56 pm** (−86% vs default 391 pm)
- N_s=10: GA **93 pm** (−77% vs default 402 pm)
- N_s=10+grad: PSO **130 pm** (−72% vs default 463 pm)
- N_s=20: GA **287 pm** (−28% vs default 396 pm)
- p=10 wybrany w 100% iteracji, spectral_spread > 0 w 87.5%

---

## 2025-12-15 | Zadanie 4 - skrypt pełnej optymalizacji

### Nowe pliki
- `scripts/Zadanie4_FullOptimization.m` - kompletny skrypt: 4 warianty × (GA + PSO), wykresy porównawcze, ~10-15 min

### Warianty
- V1: N_s=5, grad=0nm
- V2: N_s=10, grad=0nm
- V3: N_s=10, grad=0.2nm (gradient temperaturowy)
- V4: N_s=20, grad=0nm

---

## 2025-12-12 | Zadanie 4 (WP3.2+3.3) - optymalizacja GA vs PSO (5 zmiennych)

### Nowe pliki
- `src/system/optimObjective.m` - funkcja celu: MAE(deltaneff, p, Nb, n_lambdas, spectral_spread)
- `results/WP3/optimization/optimization_ga_pso.mat` - wyniki GA vs PSO

### Wyniki
5 zmiennych, N_s=10, grad=0. GA(10×8) vs PSO(10×8), oba ~25s:
- **Default** (p=8, dn=0.3e-4, Nb=8, 16λ, spread=0): MAE=381.5 pm
- **GA**: p=10, dn=4.8e-5, Nb=8, 16λ, spread=0.44nm → MAE=293.9 pm (**−23%**)
- **PSO**: p=10, dn=1.0e-5, Nb=16, 16λ, spread=0.18nm → MAE=149.2 pm (**−61%**)

PSO 2× lepszy od GA. Oba wybrały p=10 i spectral_spread > 0.

### Przyspieszenie symulacji
- `lowpass()` → FIR(64) w `runSimulation.m`: oszczędność ~0.25s/eval
- `SILENT_MODE` w generateValues/addNoise: brak niepotrzebnych printów
- Wynik: ~0.2s/eval (N_s=5) do ~0.5s/eval (N_s=10) - umożliwia optymalizację w <1 min

---

## 2025-12-10 | Zadanie 4 (WP3.2+3.3) - optymalizacja p × deltaneff (grid search)

### Nowe pliki
- `scripts/Zadanie4_Optimization.m` - grid search p×deltaneff, 2 warianty
- `results/WP3/optimization/optimization_heatmap.fig/.png` - heatmapy MAE
- `results/WP3/optimization/optimization_results.mat` - dane

### Wyniki
Grid search: p ∈ {4,6,8,10} × deltaneff ∈ [0.1e-4, 1e-4], 8 punktów. 32 symulacje w ~14s.

| Wariant | p_opt | deltaneff_opt | R_peak | MAE_opt |
|---------|-------|---------------|--------|---------|
| N_s=10, grad=0nm | 10 | 3.57e-5 | 4.7% | **269 pm** |
| N_s=10, grad=0.2nm | 10 | 2.29e-5 | 1.9% | **94 pm** |

Wniosek: p=10 zawsze optymalne. Z gradientem optimum przesuwa się na niższe deltaneff (mniej crosstalku = lepsza detekcja przesunięć). Poprawa vs default (p=8, dn=0.3e-4): V1: 381→269 pm (29%), V2: 437→94 pm (78%).

---

## 2025-11-28 | Analiza PIN vs APD - praktyczne wyniki (WP3)

### Nowe/zaktualizowane pliki
- `scripts/Zadanie_PINvsAPD.m` - kompletny skrypt: sweep gain + sweep moc optyczna
- `results/WP3/sensitivity/pin_vs_apd.fig/.png` - MAE vs gain
- `results/WP3/sensitivity/pin_vs_apd_power.fig/.png` - **MAE vs moc optyczna (Thorlabs)**
- `results/WP3/sensitivity/pin_vs_apd_power_results.mat` - dane

### Wyniki kluczowe
Parametry Thorlabs: PIN (R=0.9 A/W, NEP=15 pW), APD (M=10, NEP=2 pW input-referred).

**Punkt przecięcia ≈ −20 dBm (10 µW)**:
- P < −20 dBm: APD lepszy (do 80 pm różnicy przy −30 dBm / 1 µW)
- P > −20 dBm: PIN lepszy lub równy
- P > −15 dBm: oba limitowane przez sidelobes ≈ 393 pm - detektor nie ma znaczenia

Model excess noise: F(M) = k·M + (1−k)·(2−1/M), k=0.7 (InGaAs).

---

## 2025-11-20 | Migracja wszystkich skryptów na runSimulation

### Zmodyfikowane pliki
- `scripts/Zadanie1_MultiBounce.m` - przepisany na `runSimulation()` + `defaultParams()`
- `scripts/Zadanie2_TemperatureShift.m` - przepisany na `runSimulation()` + `defaultParams()`
- `scripts/Zadanie3_SensitivityAnalysis.m` - przepisany na `runSimulation()` + `defaultParams()`
- `scripts/Zadanie6_BenchmarkCodesDenoising.m` - przepisany na `runSimulation()` + `defaultParams()`
- `scripts/WP2_SystemSimulation_v2.m` - przepisany na `runSimulation()` + `defaultParams()`

### Efekt
Wszystkie skrypty korzystają z jednej funkcji symulacyjnej (`runSimulation`) i jednego zestawu parametrów domyślnych (`defaultParams`). Czas benchmarku 6 konfiguracji: 5.1s (vs ~120s wcześniej).

---

## 2025-11-15 | Refaktoryzacja: runSimulation + defaultParams

### Nowe pliki
- `src/system/runSimulation.m` - czysta funkcja symulacyjna, 1 struct input/output, opcjonalne ploty (`params.show_plots`, domyślnie false), automatyczny oblicz MAE/MAX. Czas: ~2.4s (vs ~20s z plotami).
- `src/system/defaultParams.m` - domyślne parametry systemu w jednym miejscu.

### Zmodyfikowane pliki
- `src/system/simulationFunction.m` - dodana flaga `fbg.silent` do wyłączania plotów
- `src/plots/plotXcovs.m` - dodane `'theme', 'light'`
- `src/plots/plotGratingsSnr.m` - dodane `'theme', 'light'`

---

## 2025-10-25 | Zadanie 3 (WP3.1) - analiza wrażliwości

### Nowe pliki
- `scripts/Zadanie3_SensitivityAnalysis.m` - skrypt z 5 przemiataniami parametrów
- `results/WP3/sensitivity/sensitivity_all.fig` + `.png` - wykres zbiorczy (6 paneli)
- `results/WP3/sensitivity/sensitivity_results.mat` - dane liczbowe

### Wyniki
Ranking parametrów wg wpływu na MAE: **p** (ΔMAE=411 pm) >> gradient T (206 pm) >> Δn_eff (47 pm) ≈ N_s (41 pm) ≈ NEP (40 pm). Kluczowe: p=10 daje MAE=73 pm (5× poprawa vs p=8). NEP ma wpływ tylko przy słabym sygnale (R<1%, N_s=10).

---

## 2025-09-15 | Zadanie 6 (WP2.2 + WP2.3) - analiza porównawcza kodów i odszumiania

### Nowe pliki
- `scripts/Zadanie6_BenchmarkCodesDenoising.m` - skrypt analizy porównawczej
- `results/WP2/benchmark/benchmark_results.mat` - wyniki
- `results/WP2/benchmark/benchmark_results.fig/.png` - wykresy porównawcze
- `results/WP2/benchmark/noise_vs_sidelobes.fig/.png` - dowód: floor = sidelobes
- `results/WP2/benchmark/filtering_vs_code_length.fig/.png` - filtracja vs długość kodu
- `results/WP2/benchmark/code_length_impact.fig/.png` - MAE vs p
- `results/WP2/benchmark/denoising_comparison.fig/.png` - odszumianie przed/po

### Wyniki benchmarku (5 FBG, p=8, deltaneff=0.3e-4)
- **Kasami > Gold**: MAE 381 vs 487 pm (różnica 27%)
- **Odszumianie nieskuteczne**: sgolay/wavelet/Gauss/RRC/TVD - żadne nie poprawia MAE
- **Odkrycie**: 100% floor korelacji to boczne listki kodowe (szum pomijalny: ΔSNR = −0.16%)
- **Długość kodu**: p=10 (L=1023) daje MAE=73 pm - jedyny skuteczny parametr poprawy

---

## 2025-07-10 | Zadanie 2 (WP2) - przesunięcie temperaturowe siatek

### Nowe pliki
- `scripts/Zadanie2_TemperatureShift.m` - skrypt demonstracyjny

### Zmodyfikowane pliki
- `src/fbg/generateAllGratings.m` - nowy parametr `fbg.lambda_shifts` (wektor przesunięć λ_B [m] per siatka). Domyślnie `zeros(1, N_s)`. Realizacja: modyfikacja periodu siatki. Gdy shifts ≠ 0, optymalizacja replikacji widm wyłączona.

### Walidacja
- Gradient 0→200 pm: MAE=436 pm (vs 382 pm baseline, +14%)
- Gradient 0→500 pm: MAE=587 pm (+54%)
- Gradient 0→2 nm: siatki poza zakresem lasera - fizycznie poprawne

---

## 2025-05-20 | Zadanie 1 (WP2) - model wielokrotnych odbić + refaktoryzacja

### Nowe pliki
- `src/fbg/generateBraggGratingTanh.m` - analityczny model FBG (teoria modów sprzężonych, siatka jednorodna). Przyspieszenie generacji widm >6000× vs TMM.
- `scripts/WP2_SystemSimulation_v2.m` - przebudowany skrypt symulacyjny (pipeline rozwinięty, zmienne w workspace)
- `scripts/Zadanie1_MultiBounce.m` - skrypt demonstracyjny

### Zmodyfikowane pliki
- `src/fbg/addAllGratings.m` - **model wielokrotnych odbić**:
  - `fbg.max_bounces = 1`: primary only (domyślne, kompatybilne wstecz)
  - `fbg.max_bounces = 3`: dodaje 3-bounce ghosts na pozycjach D(j)−D(k)+D(m)
  - Usunięta stara funkcja `crossOld()` (zastąpiona fizycznym modelem)
- `src/fbg/generateAllGratings.m` - wybór modelu FBG (`'tanh'`/`'tmm'`), `parfor` → `for`
- `src/fbg/generateBraggGrating.m` - usunięta zależność od Fuzzy Logic Toolbox
- `src/system/generateValues.m` - `L_max` uwzględnia pozycje ghostów (parametr `max_bounces`)
- `src/system/simulationFunction.m` - przekazuje `max_bounces` do `generateValues`

### Walidacja
- `max_bounces=1`: identyczne wyniki jak baseline (MAE=388.98 pm)
- `max_bounces=3`, R≈30%: ghost/primary=10%, ghosty na pozycjach FBG3/4/5 + 700–1000m
