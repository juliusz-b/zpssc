[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15089768.svg)](https://doi.org/10.5281/zenodo.15089768)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
![GitHub License](https://img.shields.io/github/license/juliusz-b/zpssc)


# 🔬 Zaawansowane techniki przetwarzania sygnałów w światłowodowych sieciach czujnikowych (ZPSSC)
 
## 📋 Opis projektu

W niniejszym repozytorium znajduje się implementacja symulatora do światłowodowych systemów czujnikowych bazujących na multipleksacji kodowej wraz z elementami optymalizacji reflektancji siatek.

> **Uwaga:** Repozytorium jest ciągle aktualizowane. Część funkcji jest testowa i w obecnej formie nie posiada szczegółowych opisów lub konwencja składniowa jest przemieszana.

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
│   ├── signals/           # Funkcje związane z przetwarzaniem sygnałów
│   └── system/            # Funkcje dotyczące symulacji systemu czujnikowego
├── scripts/               # Skrypty w MATLAB
└── tests/                 # Testowe funkcje i/lub skrypty
```

## 🚀 Jak uruchomić kody

### Wymagania

- MATLAB 2022b lub nowszy
- Signal Processing Toolbox
- Communications Toolbox

### Uruchomienie

```matlab
% Add all directories to path
AddAllSubfolders;

% Run a basic simulation
run('scripts/WP2_CodeAnalysis.m');
```

## 📊 Etapy projektu

Projekt podzielono na następujące etapy:

### WP1: Analiza i przegląd literatury ✅

#### Zakończone działania:
- Przeprowadzono dogłębną analizę rodzin kodowych: Kasamiego, PRBS, Randi, Golda, OOC, Sidelnikova, pary Golaya i sekwencji chaotycznych
- Stworzono zestaw skryptów do testowania różnych scenariuszy symulacyjnych
- Określono wstępne parametry pracy systemu kodowego (minimalne pasmo odbiornika: 20 MHz)
- Wstępne analizy wykazały lepszą detekowalność sekwencji Kasamiego w porównaniu do pozostałych sekwencji

#### Wnioski:
- Filtracja dolnoprzepustowa powoduje zmianę PSNR
- Konieczne jest korzystanie z wolniejszych kodów przy detektorze o stosunkowo małym paśmie
- Parametry pasma odbiornika są ściśle zależne od: odległości między czujnikami, wymaganej szybkości ściągania danych i długości stosowanych kodów


### WP2: Stworzenie symulatora systemu czujnikowego wraz z makietą pomiarową 🔄

#### W trakcie realizacji:
- Utworzono podstawowy symulator systemu czujnikowego
- Zaimplementowano model sieci z możliwością definiowania punktów pomiarowych
- Przygotowano środowisko do analizy metod detekcji czujników z wykorzystaniem korelacji
- Trwają prace nad implementacją metod redukcji szumów i poprawy SNR

#### Planowane działania:
- Wykonanie stanowiska laboratoryjnego do pomiaru sieci czujnikowych z możliwością zmiany temperatury wybranego czujnika światłowodowego
- Implementacja zaawansowanych metod detekcji, w tym SIC
- Analiza metod redukcji szumów, takich jak: TVD (Total Variation Denoising), adaptacyjne filtrowanie, filtracja pasmowo-przepustowa i filtr Savitzkiego-Golaya
### WP3: Optymalizacja projektowania systemów czujnikowych 📝

#### Planowane działania:
- Określenie kluczowych parametrów systemu podlegających optymalizacji
- Wybór kryteriów optymalizacyjnych i funkcji celu
- Analiza i wybór algorytmów optymalizacyjnych
- Przeprowadzenie symulacji porównawczych "przed" i "po" optymalizacji
- Stworzenie dedykowanego modelu optymalizacyjnego dla systemów z multipleksacją kodową

### WP4: Wykonanie testów laboratoryjnych 📝

#### Planowane działania:
- Zaprojektowanie układu wzmacniacza do fotodetektora o wysokim wzmocnieniu i niskich szumach własnych
- Przygotowanie stanowiska pomiarowego
- Wykonanie i walidacja modelu optymalizacyjnego
- Testy praktyczne z wykorzystaniem siatek Bragga o zoptymalizowanych parametrach

## 📜 License

Sprawdź plik LICENSE.


## 🤝 Współpraca i kontakt

Jesteśmy otwarci na współpracę z jednostkami naukowymi i podmiotami przemysłowymi zainteresowanymi zastosowaniem światłowodowych sieci czujnikowych.

Jeśli jesteś zainteresowany/-a:
- implementacją systemu czujnikowego w swojej aplikacji
- współpracą badawczą w dziedzinie fotoniki i systemów czujnikowych
- wymianą doświadczeń w zakresie przetwarzania sygnałów

Skontaktuj się z nami poprzez:
- Email: [juliusz.bojarczuk@pw.edu.pl](mailto:juliusz.bojarczuk@pw.edu.pl)
- GitHub: Otwórz Issue lub Pull Request w tym repozytorium

## 🙏 Acknowledgments

Projekt finansowany z programu "Perły Nauki" Ministra Nauki i Szkolnictwa Wyższego (nr umowy: PN/01/0321/2022).

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

