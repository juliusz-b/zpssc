<p align="center">
  <img src="docs/perly_nauki_logo.png" alt="Perły Nauki - Ministerstwo Nauki i Szkolnictwa Wyższego" width="300">
</p>

<p align="center">
  <strong>Projekt finansowany w ramach programu Perły Nauki MEiN</strong><br>
  Umowa nr PN/01/0321/2022
</p>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15089768.svg)](https://doi.org/10.5281/zenodo.15089768)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2025b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
![GitHub License](https://img.shields.io/github/license/juliusz-b/zpssc)
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=juliusz-b/zpssc)


# Zaawansowane techniki przetwarzania sygnałów w światłowodowych sieciach czujnikowych (ZPSSC)
 
## 📋 Opis projektu

Repozytorium zawiera symulator światłowodowych sieci czujnikowych z multipleksacją kodową (CDM) i optymalizacją siatek Bragga.

> **Uwaga:** Część funkcji jest w fazie testowej i może nie posiadać pełnej dokumentacji.

## 🎯 Tematyka projektu

Główne założenia projektu:

1. Analiza kodów możliwych do zastosowania w systemach czujnikowych
2. Opracowanie symulatora do sieci czujnikowych
3. Optymalizacja doboru siatek Bragga do systemu
4. Opracowanie elektronicznych układów wykorzystywanych w analizowanych systemach czujnikowych

### Schemat koncepcyjny systemu:

```
                  ┌─────────────────┐
                  │                 │
┌───────────┐     │   Światłowód    │    ┌───────────┐    ┌───────────────┐
│ Generator │     │  z czujnikami   │    │ Fotodetek-│    │ Zaawansowane  │
│  kodów    ├────►│   (FBG1...n)    ├───►│    tor    ├───►│ przetwarzanie │
│ optycznych│     │                 │    │           │    │   sygnałów    │
└───────────┘     │                 │    └───────────┘    └───────────────┘
                  └─────────────────┘
```


## 📁 Struktura repozytorium

```
.
├── src/                   # Kody źródłowe
│   ├── codes/             # Funkcje do generowania sekwencji kodowych
│   ├── fbg/               # Funkcje związane z symulacją siatek Bragga
│   ├── opt_source/        # Funkcje związane z symulacją lasera VCSEL
│   ├── plots/             # Funkcje do wyświetlania wyników
│   ├── signal/            # Funkcje związane z przetwarzaniem sygnałów
│   └── system/            # Funkcje dotyczące symulacji systemu czujnikowego
├── scripts/               # Skrypty w MATLAB
└── tests/                 # Testowe funkcje i/lub skrypty
```

## 🚀 Jak uruchomić kody

### Wymagania

- MATLAB R2025b lub nowszy
- Signal Processing Toolbox
- Communications Toolbox
- Global Optimization Toolbox
- DSP System Toolbox

### Uruchomienie

```matlab
% Add all directories to path
AddAllSubfolders;

% Run a basic simulation with default parameters
params = defaultParams();
params.show_plots = true;
out = runSimulation(params);
fprintf('MAE = %.1f pm\n', out.MAE);

% Run with Gaussian fit (better accuracy)
out_gauss = runSimulationGauss(params);
fprintf('MAE (Gaussian) = %.1f pm\n', out_gauss.MAE);
```

### Skrypty demonstracyjne

| Skrypt | Opis |
|--------|------|
| `scripts/WP2_SystemSimulation_v2.m` | Podstawowa symulacja z wykresami |
| `scripts/Zadanie1_MultiBounce.m` | Model wielokrotnych odbić |
| `scripts/Zadanie2_TemperatureShift.m` | Wpływ gradientu temperaturowego |
| `scripts/Zadanie3_SensitivityAnalysis.m` | Analiza wrażliwości (5 parametrów) |
| `scripts/Zadanie4_FullOptimization.m` | Optymalizacja GA vs PSO |
| `scripts/Zadanie4_GaussOptimization.m` | Optymalizacja z Gaussian fit |
| `scripts/Zadanie4_Optimization.m` | Grid search p x deltaneff |
| `scripts/Zadanie5_Comparison.m` | Porównanie przed/po optymalizacji |
| `scripts/Zadanie6_BenchmarkCodesDenoising.m` | Porównanie kodów i metod odszumiania |
| `scripts/Zadanie_DFE_test.m` | Test equalizera DFE |
| `scripts/Zadanie_PINvsAPD.m` | Porównanie PIN vs APD |
| `scripts/Zadanie_TDMvsCDM_Ghosty.m` | Porównanie TDM vs CDM |

## 📊 Etapy projektu

Projekt podzielono na następujące etapy:

### WP1: Analiza i przegląd literatury ✅

#### Zakończone działania:
- Przeprowadzono szczegółową analizę rodzin kodowych: Kasamiego, PRBS, pseudolosowych, Golda, OOC, Sidelnikowa, pary Golaya i sekwencji chaotycznych
- Stworzono zestaw skryptów do testowania różnych scenariuszy symulacyjnych
- Określono wstępne parametry pracy systemu kodowego (minimalne pasmo odbiornika: 20 MHz)
- Wstępne analizy wykazały lepszą detekowalność sekwencji Kasamiego w porównaniu do pozostałych sekwencji

#### Wnioski:
- Filtracja dolnoprzepustowa wpływa na stosunek sygnału do szumu (PSNR)
- Konieczne jest korzystanie z wolniejszych kodów przy detektorze o stosunkowo małym paśmie
- Parametry pasma odbiornika są ściśle zależne od: odległości między czujnikami, wymaganej szybkości ściągania danych i długości stosowanych kodów


### WP2: Stworzenie symulatora systemu czujnikowego wraz z makietą pomiarową ✅

#### Zakończone działania:
- Symulator sieci czujnikowej z pełnym pipeline: generacja kodów, widma FBG, odbicia, szum, korelacja, detekcja
- Model wielokrotnych odbić z parametrem `max_bounces` (sygnały odbić wielokrotnych 3. rzędu)
- Model przesunięcia temperaturowego siatek (`lambda_shifts`)
- Analityczny model FBG (tanh) - przyspieszenie generacji widm o >6000x
- Porównanie kodów: Kasami lepszy od Gold o 27%
- Analiza szumu: 100% szumu korelacyjnego to deterministyczne listki boczne (nie szum detektora)
- Systematyczne porównanie 8 metod odszumiania - żadna nie dała istotnej poprawy w systemie CDM
- Porównanie PIN vs APD z parametrami Thorlabs: punkt przecięcia przy −20 dBm
- Porównanie TDM vs CDM: CDM zapewnia naturalną ochronę przed sygnałami odbić wielokrotnych
- Dopasowanie krzywej Gaussa do estymacji długości fali Bragga: poprawa 63–88% vs średnia ważona

#### Kluczowe wyniki:
- MAE (średni błąd bezwzględny) detekcji: od 391 pm (parametry domyślne) do 21,6 pm (po optymalizacji z dopasowaniem Gaussa)
- Dłuższy kod (p=10) daje 5,3× poprawę MAE
- CDM tłumi odbicia wielokrotne: degradacja <0,5% w CDM vs >10% w TDM

### WP3: Optymalizacja projektowania systemów czujnikowych ✅

#### Zakończone działania:
- Identyfikacja 5 kluczowych parametrów: p, Δn_eff, N_s, NEP (ekwiwalentna moc szumu), gradient temperaturowy
- Ranking wrażliwości: p ≫ gradient ≫ pozostałe
- Optymalizacja 5-parametrowa (GA + PSO): Δn_eff, p, N_b, liczba długości fal, rozłożenie spektralne siatek
- 4 warianty (N_s=5/10/20, z/bez gradientu), poprawa do 94,5%
- Rozłożenie spektralne siatek jako nowy parametr projektowy (wybrany w 87,5% iteracji optymalizacji)
- Porównanie przed/po optymalizacji z wykresami per wariant

#### Kluczowe wyniki:
- Optymalizacja ze średnią ważoną: N_s=5 → 56 pm, N_s=10 → 93 pm
- Optymalizacja z dopasowaniem Gaussa: N_s=5 → **21,6 pm**, N_s=10 → **49,7 pm**
- p=10 i n_lambdas=16 optymalne we wszystkich wariantach

### WP4: Wykonanie testów laboratoryjnych 📝

#### Planowane działania:
- Zaprojektowanie układu wzmacniacza do fotodetektora o wysokim wzmocnieniu i niskich szumach własnych
- Przygotowanie stanowiska pomiarowego
- Wykonanie i walidacja modelu optymalizacyjnego
- Testy praktyczne z wykorzystaniem siatek Bragga o zoptymalizowanych parametrach

## 📜 Licencja

Sprawdź plik LICENSE.


## 🤝 Współpraca i kontakt

Jesteśmy otwarci na współpracę z jednostkami naukowymi i podmiotami przemysłowymi zainteresowanymi zastosowaniem światłowodowych sieci czujnikowych.

Jeśli jesteś zainteresowany/-a:
- implementacją systemu czujnikowego w swojej aplikacji
- współpracą badawczą w dziedzinie fotoniki i systemów czujnikowych
- wymianą doświadczeń w zakresie przetwarzania sygnałów

Skontaktuj się ze mną poprzez:
- Email: [juliusz.bojarczuk@pw.edu.pl](mailto:juliusz.bojarczuk@pw.edu.pl)
- GitHub: Otwórz Issue lub Pull Request w tym repozytorium

## 🙏 Acknowledgments

Projekt finansowany z programu "Perły Nauki" Ministerstwa Nauki i Szkolnictwa Wyższego (nr umowy: PN/01/0321/2022).

## 📚 Bibliografia

[1] H.-E. Joe, H. Yun, S.-H. Jo, M. B. G. Jun, and B.-K. Min, "A review on optical fiber sensors for environmental monitoring," *Int. J. Precis. Eng. Manuf. Technol.*, vol. 5, pp. 173–191, 2018.

[2] Y. Rao, "Recent progress in applications of in-fibre Bragg grating sensors," *Opt. Lasers Eng.*, vol. 31, pp. 297–324, 1999.

[3] J. Leng and W. Ecke, "Opportunities of fiber optic sensors and their applications," *Opt. Lasers Eng.*, vol. 47, p. 1017, 2009.

[4] F. Taffoni, D. Formica, P. Saccomandi, G. Pino, and E. Schena, "Optical Fiber-Based MR-Compatible Sensors for Medical Applications: An Overview," *Sensors*, vol. 13, pp. 14105–14120, 2013.

[5] S. Silvestri and E. Schena, "Optical-Fiber Measurement Systems for Medical Applications," in *Optoelectronics - Devices and Applications*, InTech, 2011.

[6] A. D. Kersey et al., "Fiber grating sensors," *J. Lightwave Technol.*, vol. 15, pp. 1442–1463, 1997.

[7] J. Chen, B. Liu, and H. Zhang, "Review of fiber Bragg grating sensor technology," *Front. Optoelectron. China*, vol. 4, pp. 204–212, 2011.

[8] K. P. Koo, A. B. Tveten, and S. T. Vohra, "DWDM of fiber Bragg grating sensors without sensor spectral dynamic range limitation using CDMA," in *OFC/IOOC 1999 - Optical Fiber Communication Conference*, vol. 4, pp. 168–170, 1999.

[9] H. Jiang et al., "Wavelength detection of model-sharing fiber Bragg grating sensor networks using long short-term memory neural network," *Opt. Express*, vol. 27, p. 20583, 2019.

[10] C. Z. Shi, C. C. Chan, W. Jin, Y. B. Liao, Y. Zhou, and M. S. Demokan, "Improving the performance of a FBG sensor network using a genetic algorithm," *Sensors Actuators A Phys.*, vol. 107, pp. 57–61, 2003.

[11] H. Jiang, J. Chen, and T. Liu, "Multi-objective design of an FBG sensor network using an improved Strength Pareto Evolutionary Algorithm," *Sensors Actuators A Phys.*, vol. 220, pp. 230–236, 2014.

[12] A. Triana, D. Pastor, and M. Varón, "A Code Division Design Strategy for Multiplexing Fiber Bragg Grating Sensing Networks," *Sensors*, vol. 17, p. 2508, 2017.

[13] M. Gotten, S. Lochmann, A. Ahrens, E. Lindner, and J. Van Roosbroeck, "2000 Serial FBG Sensors Interrogated with a Hybrid CDM-WDM Scheme," *J. Lightwave Technol.*, vol. 38, pp. 2493–2503, 2020.

[14] B. Crockett, L. Romero Cortés, R. Maram, and J. Azaña, "Optical signal denoising through temporal passive amplification," *Optica*, vol. 9, p. 130, 2022.

[15] D. Tosi, "Review and analysis of peak tracking techniques for fiber bragg grating sensors," *Sensors*, vol. 17, art. no. 2368, 2017.

[16] Z. Zhou et al., "Optical fiber Bragg grating sensor assembly for 3D strain monitoring and its case study in highway pavement," *Mech. Syst. Signal Process.*, vol. 28, pp. 36–49, 2012.

[17] A. Barrias, J. Casas, and S. Villalba, "A Review of Distributed Optical Fiber Sensors for Civil Engineering Applications," *Sensors*, vol. 16, p. 748, 2016.

[18] D. A. Krohn, T. W. MacDougall, and A. Mendez, *Fiber Optic Sensors: Fundamentals and Applications*. Society of Photo-Optical Instrumentation Engineers (SPIE), 2014.

