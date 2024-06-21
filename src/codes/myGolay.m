function seq = myGolay(len)

[a,b] = generateGolayPair(len);
seq = [a;b];

end


function [a,b] = generateGolayPair(length_type)
% GENERATEGOLAYPAIR Generuje parę sekwencji komplementarnych Golaya
%
% Składnia:
%   [a, b] = generateGolayPair(length_type)
%
% Wejście:
%   length_type - określenie długości sekwencji:
%       - liczba: 2, 4, 8, 16, 32, 64, ... (potęgi 2)
%       - 'seed': zwraca parę bazową [+1, +1], [+1, -1]
%       - 10, 26: specjalne długości (nie będące potęgami 2)
%
% Wyjście:
%   a - pierwsza sekwencja komplementarna Golaya
%   b - druga sekwencja komplementarna Golaya
%
% Przykłady użycia:
%   [a, b] = generateGolayPair(16);    % Sekwencje o długości 16
%   [a, b] = generateGolayPair('seed'); % Bazowa para inicjująca

% Sprawdzenie argumentu wejściowego
if ischar(length_type) && strcmpi(length_type, 'seed')
    % Zwróć parę inicjującą (seed)
    a = [1, 1];
    b = [1, -1];
    return;
elseif length_type == 10
    % Specjalna para długości 10
    a = [1, 1, 1, -1, 1, 1, -1, -1, 1, -1];
    b = [1, 1, 1, 1, -1, -1, -1, 1, -1, -1];
    return;
elseif length_type == 26
    % Specjalna para długości 26
    a = [1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, -1, -1, -1, -1, 1, -1, 1, -1, 1, 1];
    b = [1, 1, 1, 1, 1, 1, -1, -1, -1, -1, 1, -1, 1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, 1, -1, -1];
    return;
elseif ~isPowerOfTwo(length_type)
    error('Długość musi być potęgą 2 (lub specjalnie obsługiwaną długością 10 lub 26)');
end

% Rekurencyjna konstrukcja sekwencji Golaya
if length_type == 2
    % Bazowa para inicjująca
    a = [1, 1];
    b = [1, -1];
else
    % Rekurencyjne wywołanie dla mniejszej długości
    [a_prev, b_prev] = generateGolayPair(length_type/2);
    
    % Konstruujemy dłuższe sekwencje zgodnie z rekurencją Golaya
    a = [a_prev, b_prev];
    b = [a_prev, -b_prev];
end
end

function result = isPowerOfTwo(n)
% Sprawdzenie, czy liczba jest potęgą 2
if n <= 0
    result = false;
else
    result = bitand(n, (n - 1)) == 0;
end
end

function [W, H] = generateHadamardMatrix(n)
% Generowanie macierzy Hadamarda o wymiarze n
if n == 1
    W = 1;
    H = 1;
    return;
end

if ~isPowerOfTwo(n)
    error('Wymiar macierzy Hadamarda musi być potęgą 2');
end

% Rekurencyjna konstrukcja macierzy Hadamarda
[W_prev, H_prev] = generateHadamardMatrix(n/2);

% Macierz Walsha (uporządkowana)
W = [W_prev, W_prev; W_prev, -W_prev];

% Macierz Hadamarda (naturalna)
H = [H_prev, H_prev; H_prev, -H_prev];
end

function [a, b] = generateGolayPairHadamard(n)
% Alternatywna metoda generowania sekwencji Golaya za pomocą macierzy Hadamarda
% Działa tylko dla długości n = 2^m

if ~isPowerOfTwo(n)
    error('Długość musi być potęgą 2');
end

% Generujemy macierz Hadamarda
[~, H] = generateHadamardMatrix(n);

% Wybieramy dwa wiersze, które tworzą parę Golaya
a = H(1, :);
b = H(2, :);
end

function golay_set = generateGolayComplementarySet(m, q)
% GENERATEGOLAYCOMPLEMENTARYSETS Generuje zbiór komplementarnych sekwencji Golaya
%
% Składnia:
%   golay_set = generateGolayComplementarySet(m, q)
%
% Wejście:
%   m - liczba sekwencji w zbiorze (zwykle potęga 2)
%   q - długość każdej sekwencji (zwykle potęga 2)
%
% Wyjście:
%   golay_set - macierz m×q zawierająca m sekwencji o długości q
%
% Uwaga:
%   Zbiory komplementarne Golaya to uogólnienie par Golaya,
%   gdzie suma autokorelacji wszystkich sekwencji w zbiorze
%   ma zerowe boczne listki.

% Zastosujemy metodę konstrukcji opartą na macierzach Hadamarda
if ~isPowerOfTwo(m) || ~isPowerOfTwo(q)
    error('Zarówno m (liczba sekwencji) jak i q (długość) powinny być potęgami 2');
end

% Generujemy macierz Hadamarda o wymiarze m
[~, H] = generateHadamardMatrix(m);

% Generujemy macierz Hadamarda o wymiarze q
[~, H_q] = generateHadamardMatrix(q);

% Inicjalizacja zbioru
golay_set = zeros(m, q);

% Wypełniamy zbiór
for i = 1:m
    golay_set(i, :) = H(i, :);
end

% Dla większych zbiorów można zastosować bardziej zaawansowane metody konstrukcji
end

function [R_a, R_b, R_sum] = computeGolayCrossCorrelation(a, b)
% Oblicza autokorelację sekwencji a i b oraz ich sumę
R_a = xcorr(a);
R_b = xcorr(b);
R_sum = R_a + R_b;
end

function [correlation_sums] = computeSetCorrelation(golay_set)
% Oblicza sumę autokorelacji dla zbioru sekwencji komplementarnych Golaya
num_sequences = size(golay_set, 1);
correlation_sum = zeros(1, 2*size(golay_set, 2)-1);

for i = 1:num_sequences
    correlation_sum = correlation_sum + xcorr(golay_set(i, :));
end

correlation_sums = correlation_sum;
end

function demoGolaySequences()
% DEMOGOLAYSEQUENCES Demonstracja właściwości sekwencji komplementarnych Golaya

% 1. Generowanie i analiza klasycznych par Golaya
lengths = [2, 4, 8, 16, 32];
figure;
sgtitle('Właściwości sekwencji komplementarnych Golaya');

for i = 1:length(lengths)
    [a, b] = generateGolayPair(lengths(i));
    [R_a, R_b, R_sum] = computeGolayCrossCorrelation(a, b);
    
    % Normalizacja dla lepszej wizualizacji
    R_a = R_a / lengths(i);
    R_b = R_b / lengths(i);
    R_sum = R_sum / (2*lengths(i));
    
    subplot(length(lengths), 3, (i-1)*3+1);
    stem(a, 'filled');
    title(['Sekwencja a, długość ', num2str(lengths(i))]);
    ylim([-1.2, 1.2]);
    
    subplot(length(lengths), 3, (i-1)*3+2);
    stem(b, 'filled');
    title(['Sekwencja b, długość ', num2str(lengths(i))]);
    ylim([-1.2, 1.2]);
    
    subplot(length(lengths), 3, (i-1)*3+3);
    lags = -(lengths(i)-1):(lengths(i)-1);
    plot(lags, R_a, 'r-', lags, R_b, 'g-', lags, R_sum, 'b-', 'LineWidth', 1.5);
    title('Autokorelacje i ich suma');
    legend('R_a', 'R_b', 'R_{suma}');
    grid on;
end

% 2. Właściwości widmowe sekwencji Golaya
n = 64;
[a, b] = generateGolayPair(n);

figure;
sgtitle('Właściwości widmowe sekwencji Golaya');

% Transformaty Fouriera
A = fft(a);
B = fft(b);

% Płaskie widmo mocy sumy kwadratów
power_sum = abs(A).^2 + abs(B).^2;

subplot(2, 2, 1);
plot(abs(A));
title('Amplituda transformaty Fouriera sekwencji a');
xlabel('Indeks częstotliwości');
ylabel('|A(k)|');
grid on;

subplot(2, 2, 2);
plot(abs(B));
title('Amplituda transformaty Fouriera sekwencji b');
xlabel('Indeks częstotliwości');
ylabel('|B(k)|');
grid on;

subplot(2, 2, 3);
plot(power_sum);
title('Suma kwadratów amplitud |A(k)|^2 + |B(k)|^2');
xlabel('Indeks częstotliwości');
ylabel('Moc');
grid on;

subplot(2, 2, 4);
hist(power_sum, 20);
title('Histogram sumy kwadratów amplitud');
xlabel('Wartość');
ylabel('Liczność');
grid on;

% 3. Demonstracja zbiorów komplementarnych Golaya
figure;
sgtitle('Zbiory komplementarne Golaya (uogólnione)');

% Generujemy zbiór 4 sekwencji o długości 16
m = 4; % Liczba sekwencji
q = 16; % Długość każdej sekwencji
golay_set = generateGolayComplementarySet(m, q);

% Wyświetlamy sekwencje
for i = 1:m
    subplot(m+1, 1, i);
    stem(golay_set(i,:), 'filled');
    title(['Sekwencja ', num2str(i)]);
    ylim([-1.2, 1.2]);
    grid on;
end

% Obliczamy sumę autokorelacji
correlation_sum = computeSetCorrelation(golay_set);

% Wyświetlamy sumę autokorelacji
subplot(m+1, 1, m+1);
lags = -(q-1):(q-1);
plot(lags, correlation_sum / m, 'LineWidth', 1.5);
title(['Suma autokorelacji ', num2str(m), ' sekwencji']);
grid on;

% Wartość szczytu powinna być równa m*q
peak_value = correlation_sum(q);
fprintf('Teoretyczna wartość szczytu: %d\n', m*q);
fprintf('Rzeczywista wartość szczytu: %f\n', peak_value);
fprintf('Maksymalna wartość bocznego listka: %f\n', max(abs(correlation_sum([1:q-1, q+1:end]))));

% 4. Porównanie sekwencji Golaya o różnych długościach
specials = [10, 26];
figure;
sgtitle('Specjalne długości sekwencji Golaya');

for i = 1:length(specials)
    [a, b] = generateGolayPair(specials(i));
    [R_a, R_b, R_sum] = computeGolayCrossCorrelation(a, b);
    
    % Normalizacja dla lepszej wizualizacji
    R_a = R_a / specials(i);
    R_b = R_b / specials(i);
    R_sum = R_sum / (2*specials(i));
    
    subplot(length(specials), 3, (i-1)*3+1);
    stem(a, 'filled');
    title(['Sekwencja a, długość ', num2str(specials(i))]);
    ylim([-1.2, 1.2]);
    
    subplot(length(specials), 3, (i-1)*3+2);
    stem(b, 'filled');
    title(['Sekwencja b, długość ', num2str(specials(i))]);
    ylim([-1.2, 1.2]);
    
    subplot(length(specials), 3, (i-1)*3+3);
    lags = -(specials(i)-1):(specials(i)-1);
    plot(lags, R_a, 'r-', lags, R_b, 'g-', lags, R_sum, 'b-', 'LineWidth', 1.5);
    title('Autokorelacje i ich suma');
    legend('R_a', 'R_b', 'R_{suma}');
    grid on;
end
end

function demoMultiCarrierSystem()
% DEMOMULTICARRIERSYSTEM Demonstracja zastosowania sekwencji Golaya w systemach OFDM
% Pokazuje jak sekwencje Golaya mogą być wykorzystane do redukcji PAPR

% Parametry systemu
N_fft = 128;           % Rozmiar FFT (liczba podnośnych)
N_cp = 16;             % Długość prefiksu cyklicznego
N_sym = 10;            % Liczba symboli OFDM
modulation = 'QPSK';   % Typ modulacji

% Generowanie sekwencji Golaya o długości FFT/2
[a_golay, b_golay] = generateGolayPair(N_fft/2);

% Tworzenie mapowania QPSK
qpsk_map = [1+1i, 1-1i, -1+1i, -1-1i] / sqrt(2);

% Generowanie danych losowych
data_bits = randi([0, 3], 1, N_fft * N_sym);
data_qpsk = qpsk_map(data_bits + 1);

% Kształtowanie sygnału za pomocą sekwencji Golaya - zestaw 1
data_shaped_1 = zeros(1, N_fft * N_sym);
for i = 1:N_sym
    idx = (i-1)*N_fft + (1:N_fft);
    
    % Podziel dane na dwie części
    half_idx = N_fft/2;
    data1 = data_qpsk(idx(1:half_idx));
    data2 = data_qpsk(idx(half_idx+1:end));
    
    % Konstruowanie symbolu OFDM przy użyciu sekwencji Golaya
    golay_mapped = [data1 .* a_golay, data2 .* b_golay];
    data_shaped_1(idx) = golay_mapped;
end

% Kształtowanie sygnału bez sekwencji Golaya (referencja) - zestaw 2
data_shaped_2 = data_qpsk;

% Transformacja IFFT i dodanie prefiksu cyklicznego
ofdm_signal_1 = zeros(1, (N_fft + N_cp) * N_sym);
ofdm_signal_2 = zeros(1, (N_fft + N_cp) * N_sym);

for i = 1:N_sym
    data_idx = (i-1)*N_fft + (1:N_fft);
    ofdm_idx = (i-1)*(N_fft + N_cp) + (1:N_fft + N_cp);
    
    % IFFT dla zestawu 1 (z sekwencjami Golaya)
    ifft_data_1 = ifft(data_shaped_1(data_idx), N_fft);
    
    % IFFT dla zestawu 2 (bez sekwencji Golaya)
    ifft_data_2 = ifft(data_shaped_2(data_idx), N_fft);
    
    % Dodanie prefiksu cyklicznego
    ofdm_symbol_1 = [ifft_data_1(end-N_cp+1:end), ifft_data_1];
    ofdm_symbol_2 = [ifft_data_2(end-N_cp+1:end), ifft_data_2];
    
    % Zapisanie w sygnale OFDM
    ofdm_signal_1(ofdm_idx) = ofdm_symbol_1;
    ofdm_signal_2(ofdm_idx) = ofdm_symbol_2;
end

% Obliczanie PAPR dla obu sygnałów
papr_1 = calculatePAPR(ofdm_signal_1);
papr_2 = calculatePAPR(ofdm_signal_2);

% Wizualizacja wyników
figure;
sgtitle('Zastosowanie sekwencji Golaya w systemach OFDM do redukcji PAPR');

% Fragment sygnału OFDM
subplot(2, 1, 1);
t = 1:(N_fft + N_cp) * 2;
plot(t, abs(ofdm_signal_1(t)), 'b-', t, abs(ofdm_signal_2(t)), 'r--', 'LineWidth', 1.5);
title('Amplituda sygnału OFDM (fragment)');
legend('Z sekwencją Golaya', 'Bez sekwencji Golaya');
xlabel('Próbka');
ylabel('Amplituda');
grid on;

% Porównanie PAPR dla obu sygnałów
subplot(2, 1, 2);
[f_1, x_1] = ecdf(papr_1);
[f_2, x_2] = ecdf(papr_2);

semilogy(x_1, 1-f_1, 'b-', x_2, 1-f_2, 'r--', 'LineWidth', 1.5);
title('CCDF funkcji PAPR');
legend(['Z sekwencją Golaya (śr. PAPR = ', num2str(mean(papr_1), '%.2f'), ' dB)'], ...
       ['Bez sekwencji Golaya (śr. PAPR = ', num2str(mean(papr_2), '%.2f'), ' dB)']);
xlabel('PAPR [dB]');
ylabel('Pr(PAPR > abscysa)');
grid on;

fprintf('Średnie PAPR z sekwencją Golaya: %.2f dB\n', mean(papr_1));
fprintf('Średnie PAPR bez sekwencji Golaya: %.2f dB\n', mean(papr_2));
fprintf('Redukcja PAPR: %.2f dB\n', mean(papr_2) - mean(papr_1));
end

function papr = calculatePAPR(signal)
% Oblicza PAPR (Peak-to-Average Power Ratio) w dB
% for various windows of the signal

window_size = 256;
num_windows = floor(length(signal) / window_size);
papr = zeros(1, num_windows);

for i = 1:num_windows
    idx = (i-1)*window_size + (1:window_size);
    sig_window = signal(idx);
    
    peak_power = max(abs(sig_window).^2);
    avg_power = mean(abs(sig_window).^2);
    
    papr(i) = 10 * log10(peak_power / avg_power);
end
end

function result = isGolayPair(a, b)
% Sprawdza, czy sekwencje a i b tworzą parę Golaya
[~, ~, R_sum] = computeGolayCrossCorrelation(a, b);
center_idx = length(R_sum) / 2 + 0.5;
sidelobes = R_sum([1:center_idx-1, center_idx+1:end]);

% Tolerancja numeryczna
tol = 1e-10;
if max(abs(sidelobes)) < tol
    result = true;
    fprintf('Sekwencje tworzą parę komplementarną Golaya.\n');
else
    result = false;
    fprintf('Sekwencje NIE tworzą pary komplementarnej Golaya.\n');
    fprintf('Maksymalna wartość bocznego listka: %e\n', max(abs(sidelobes)));
end
end