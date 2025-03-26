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

Przeanalizowano kody Kasamiego, PRBS, Randi, Golda, OOC, Sidelnikova, pary Golaya, sekwencje chaotyczne. Filtracja dolnoprzepustowa powoduje zmianę PSNR. Konieczne jest więc korzystanie z wolniejszych kodów przy detektorze o stosunkowo małym paśmie.

### WP2: Stworzenie symulatora systemu czujnikowego wraz z makietą pomiarową 🔄

Utworzono prosty symulator. Trwają prace nad modyfikacją skryptów i aktualizacją funkcjonalności.

### WP3: Optymalizacja projektowania systemów czujnikowych 📝

*TODO*

### WP4: Wykonanie testów laboratoryjnych 📝

*TODO*

## 📜 License

Sprawdź plik LICENSE.

## 🙏 Acknowledgments

Projekt finansowany z programu "Perły Nauki" Ministra Nauki i Szkolnictwa Wyższego (nr umowy: PN/01/0321/2022).
