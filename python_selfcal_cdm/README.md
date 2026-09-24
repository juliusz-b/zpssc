# python_selfcal_cdm - samokalibrująca interrogacja CDM (badanie symulacyjne)

Lekkie, odtwarzalne symulacje w Pythonie do artykułu JCR (EN) o samokalibrującej
się interrogacji CDM sieci siatek Bragga z tanim, bezpośrednio
modulowanym, przestrajalnym laserem VCSEL. Uzupełnienie symulatora MATLAB (`../`). Kod i opisy po angielsku.

> **STAN:** po rundzie 6 recenzji (major revision), w rewizji. Handoff i decyzje:
> `../../Realizacja zadań/Artykul_JCR_2026/STAN_ARTYKULU_2026-06.md` oraz `RECENZJE_I_DECYZJE_2026.md`.

## Wymagania
- Python 3.10+, `numpy`, `scipy`, `matplotlib` (`pip install numpy scipy matplotlib`)

## Pliki
| Plik | Opis |
|------|------|
| `common.py` | wspólna fizyka + parametry sprzętu: kody Gold/Kasami, widma FBG (Gauss/tanh²), model chirpu i konwersji FM→AM, estymatory piku (Gauss/centroid), nieliniowa charakterystyka przestrajania VCSEL, linijka MZI |
| `s0_schematics.py` | schematy poglądowe: diagram blokowy + odwzorowanie długość fali-czas |
| `s1_codes.py` | właściwości korelacyjne kodów; podłoga listków vs długość kodu (Fig.3) |
| `s2_despread.py` | korelacyjne rozdzielenie (despreading) nakładających się siatek + kanały referencyjne (Fig.4, waterfall) |
| `s3_residual.py` | rezyduum błędu λ_B po korekcie referencjami vs chirp/FWHM i niedopasowanie nachylenia (Fig.5) |
| `s4_baselines.py` | strategie kalibracji: brak / MZI / referencje / referencje+MZI (Fig.6) |
| `s5_ghosts.py` | odporność punktu odniesienia na odbicia wielokrotne (Fig.7) |
| `s6_experiment.py` | predykcja makiety (FBGS DTG-A3A4, 1545 nm, 250 pm, R=10%) (Fig.8) |
| `s7_axis_nonlinearity.py` | kalibracja nieliniowej osi λ(V) HCG; szerokie pasmo vs wąska makieta (Fig.9) |
| `s_setup_figure.py` | schemat stanowiska eksperymentalnego (Fig.10) |
| `s8_noise.py` | podłoga szumowa: jitter centroidu vs SNR i reflektancja |
| `s9_coherent.py` | długość koherencji z linii CW; przypadek koherentny vs niekoherentny |
| `s10_control.py` | **kontrola nowości**: referencja ko-kodowana vs zwykła (decydujące) |
| `s11_sensitivity.py` | wrażliwość rezyduum na rozrzut niesymetrii i kształt (tanh² vs sinc²) |
| `s12_capacity.py` | **pojemność sieci**: błąd λ_B vs liczba siatek z rozkładem na mechanizmy (przysłanianie widmowe, echa pozorne, przesłuch kodowy, szum); rozstaw równomierny vs losowy; sekwencyjna korekta przysłaniania; pojemność vs reflektancja |
| `s13_chiprate.py` | **tor akwizycji**: rozdzielczość przestrzenna vs szybkość chipowa, tolerancja jitteru, budżet próbkowania ekwiwalentnego (ETS), rozdzielczość przetwornika |
| `s14_code_families.py` | **rodziny kodów**: m-sekwencje vs Gold vs Kasami vs pary Golaya przy przestrajanym źródle nadającym jeden kod naraz (liczy się autokorelacja, nie korelacja wzajemna) |
| `s15_budget.py` | **zbiorczy budżet błędu** dla dwóch konfiguracji (makieta w zakupionej konfiguracji, sieć zaprojektowana) + ile referencji się opłaca |
| `s16_principle.py` | rysunki poglądowe: zasada działania (układ + mapa opóźnienie-długość fali) i dwa mechanizmy błędu w sieci (echa pozorne, przesłuch kodowy) |
| `s17_coding.py` | co robi kod w dziedzinie czasu: przebieg, suma powrotów, wynik korelacji |
| `s18_source.py` | chirp wywołany kodem (składowa adiabatyczna + przerzut), jądro p(δ), konwersja FM-AM na zboczu, rezyduum vs liczba referencji |
| `s19_deshadow.py` | przysłanianie widmowe i rekurencyjna korekta: schemat, linia przed/po, zysk i granica |
| `s20_inversion.py` | **czego korekta musi znać**: wariant bez założenia kształtu, tolerancja skali amplitudy, rozróżnianie echa pozornego od siatki po szerokości linii, dwie siatki w jednym binie opóźnienia |
| `s21_cdmwdm.py` | **hybryda CDM-WDM**: adresowanie pasmo x bin opóźnienia, cykl pomiarowy w czasie, dowód numeryczny że pasmowanie resetuje mechanizmy błędu (32 w 1 paśmie: 60,1 pm; 4x8: 22,1; 8 w 1 paśmie: 22,6) |
| `s22_codelength.py` | **długość kodu N jako zmienna projektowa**: biny/zasięg, sufit przesłuchu 0,37N (potwierdzony dla N=63..511), czas ramki i odświeżanie, statystyka kolizji binów (paradoks dnia urodzin) i oczekiwana liczba „trudnych par" ~0,26·K²/N |
| `s23_theory.py` | **teoria w zamkniętej formie, zweryfikowana**: prawo A (bias przysłaniania pary, dokładne do 1. rzędu w R, maks. 0,81·R·sigma przy odstrojeniu sigma·sqrt(3/2)), prawo B (średni bias wielodostępu, nieparzysty, liniowy w (K-1)/N, zdejmowany przez referencje 4x), nierówność C (K·f_r <= eps·f_s/(c_L·n_s·M), N się skraca) |
| `s24_storyboard.py` | **rysunek fabularny**: życie jednego okresu kodu - światłowód z osią odległości nad pasami czasowymi w tej samej skali (t=2nz/c), echa startują pod swoimi siatkami, piki korelacji tamże |
| `test_selfcal.py` | testy: granice korelacyjne kodów, estymatory, przypadki graniczne fizyki, algebra opóźnień ech pozornych, wzory akwizycji (38 testów, `python test_selfcal.py`) |

## Uruchomienie
```bash
python3 s1_codes.py       # analogicznie s2..s15; figury -> figs/
python3 test_selfcal.py   # 30 testów, kończy się kodem != 0 przy błędzie
```

## Kluczowe wyniki
- Kody Gold N=127: korelacja wzajemna **0,134 = granica 17/127**; podłoga listków maleje z długością kodu.
- Korelacyjne rozdzielenie CDM (idealny tor): MAE ≈ **9 pm**.
- Referencje usuwają wspólny offset chirpu → **podłoga rezyduum ~kilka pm** (niedopasowanie nachylenia). Surowy błąd FM-AM **skaluje się ze chirpem** (chirp do zmierzenia na OSA; NIE zakładać konkretnej wartości).
- Oś λ(V) (s7): szerokie pasmo wymaga 3 ref lub MZI; wąska makieta: 2 ref wystarczą (~9 pm).
- Echa pozorne: punkt odniesienia (Gauss) niezmienniczy; koherencja (s9): L_coh ≈ 0,22 m ≪ 4 m → **niekoherentnie**.
- Szum (s8): przy SNR ≥16 dB systematyka rządzi.
- **Kontrola nowości (s10): ko-kodowanie referencji to zaleta SYSTEMOWA (bez dodatkowego sprzętu), NIE poprawa dokładności**, przy innych kanałach nawet szkodzi (cross-talk). Zmienia tezę pracy.

### Runda projektowa (s12-s15, sierpień 2026)
- **Przysłanianie widmowe rządzi przy silnych siatkach**: przy R=10% już 3 siatki dają 8,1 pm, a 96 siatek 179 pm. Sekwencyjna korekta przysłaniania (peeling) sprowadza makietę do 0,6 pm i podnosi pojemność z 4 do 17 siatek przy R=10%.
- **Rozstaw siatek ma znaczenie strukturalne**: przy rozstawie równomiernym KAŻDE echo pozorne trzeciego rzędu mieszczące się w rozpiętości szeregu siatek trafia w zajęty bin opóźnienia; rozstaw losowy zostawia ok. 12%. Przy K=32 to różnica 45,0 vs 10,4 pm.
- **m-sekwencje biją kody Gold** przy przestrajanym źródle (jeden kod naraz): listek autokorelacji 1/N vs 17/N, błąd 6,5 vs 28,1 pm przy K=32. Pary Golaya znoszą listek całkowicie kosztem dwukrotnego czasu akwizycji.
- **Sufit przesłuchu kodowego skaluje się jak K/N** i nie zależy od reflektancji: ok. 0,21·K pm przy N=127, czterokrotnie mniej przy N=511.
- **Rozdzielczość**: kryterium Δz = c/(2nB) potwierdzone numerycznie (4 m → 25 Mchip/s, 2 m → 50 Mchip/s). Tolerancja jitteru σ_t < 0,2 okresu chipu.
- **ETS zamyka spór o ADC**: pełny rekord przy 25 Mchip/s to 147 µs przy ADC 3,6 MS/s i 31 µs przy 20 MS/s. Kwantyzacja jest nieistotna (0,34 pm przy 8 bitach i 2% wypełnienia skali).
- **Korekta przysłaniania nie wymaga kształtu linii (s20)**: transmisję da się zbudować punkt po punkcie z samego pomiaru, `T = (1 - g·S)²`. Na siatkach tanh² wariant bezkształtowy daje 1,90 pm wobec 3,67 pm dla dopasowania gaussowskiego. Potrzebna jest tylko skala amplitudy i to z grubsza: korekta pomaga dla g od 0,5 do 1,6.
- **Granica K ≈ 24 przy R=10% to nieprzezroczystość, nie chciwy algorytm**: ostatnia siatka w szeregu 32 siatek dostaje 700 razy słabszy sygnał, a iteracja po całym szeregu zmienia wynik o mniej niż 2 pm.
- **Echo pozorne da się rozpoznać (s20)**: linia węższa o √3 (zmierzone 0,571 wobec 0,577 z teorii), amplituda ∝ R³. Jeden próg na szerokości daje 95% zbalansowanej trafności.
- **Dwie siatki w jednym binie opóźnienia (s20)**: dopasowanie dwuskładnikowe rozdziela je od pół szerokości linii (9,2 pm), przy pełnej szerokości 0,55 pm. Poniżej pół szerokości zawodzi. To pozwala wymieniać szybkość chipową na separację widmową.
- **Budżet błędu (s15)**: makieta w zakupionej konfiguracji 9,9 pm (dominuje przysłanianie 8,1 i chirp 4,8), sieć zaprojektowana 3,3 pm. Trzy referencje tną rezyduum chirpu 24-krotnie względem braku referencji, czwarta daje już niewiele.

Jednostki: częstotliwość optyczna w GHz, 1 GHz ≈ 8 pm przy 1550 nm.
