%demoChaoticSequences()



function [sequence] = myChaotic(type, length, params, initial_state, transient)
% GENERATECHAOTICSEQUENCE Generuje sekwencję chaotyczną określonego typu
%
% Składnia:
%   [sequence] = generateChaoticSequence(type, length, params, initial_state, transient)
%
% Wejście:
%   type - typ sekwencji chaotycznej: 'logistic', 'tent', 'henon', 'lorenz', 'bernoulli'
%   length - długość generowanej sekwencji
%   params - parametry dla wybranego modelu chaotycznego (struktura lub wektor)
%   initial_state - początkowy stan układu (domyślnie losowy)
%   transient - długość przebiegu przejściowego do pominięcia (domyślnie 100)
%
% Wyjście:
%   sequence - wygenerowana sekwencja chaotyczna
%
% Przykłady użycia:
%   seq = generateChaoticSequence('logistic', 1000, 3.99, 0.5, 100);
%   seq = generateChaoticSequence('lorenz', 1000, [10, 28, 8/3], [1, 1, 1], 500);

% Sprawdź liczbę argumentów wejściowych
if nargin < 2
    error('Wymagane są co najmniej argumenty: typ i długość');
end

% Domyślne wartości
if nargin < 5
    transient = 100;
end

% Inicjalizacja sekwencji
sequence = zeros(1, length);

% Generowanie sekwencji w zależności od typu
switch lower(type)
    case 'logistic'
        % Parametry dla odwzorowania logistycznego
        if nargin < 3 || isempty(params)
            r = 3.99; % Wartość domyślna dla zachowania silnego chaosu
        else
            r = params;
        end
        
        % Stan początkowy
        if nargin < 4 || isempty(initial_state)
            x = rand(); % Losowa wartość między 0 a 1
        else
            x = initial_state;
        end
        
        % Pomiń stany przejściowe
        for i = 1:transient
            x = r * x * (1 - x);
        end
        
        % Generowanie sekwencji
        for i = 1:length
            x = r * x * (1 - x);
            sequence(i) = x;
        end
        
    case 'tent'
        % Parametry dla odwzorowania namiotowego (tent map)
        if nargin < 3 || isempty(params)
            mu = 1.99; % Wartość domyślna bliska 2 dla silnego chaosu
        else
            mu = params;
        end
        
        % Stan początkowy
        if nargin < 4 || isempty(initial_state)
            x = rand(); % Losowa wartość między 0 a 1
        else
            x = initial_state;
        end
        
        % Pomiń stany przejściowe
        for i = 1:transient
            if x < 0.5
                x = mu * x;
            else
                x = mu * (1 - x);
            end
        end
        
        % Generowanie sekwencji
        for i = 1:length
            if x < 0.5
                x = mu * x;
            else
                x = mu * (1 - x);
            end
            sequence(i) = x;
        end
        
    case 'henon'
        % Parametry dla odwzorowania Hénona
        if nargin < 3 || isempty(params)
            a = 1.4;
            b = 0.3;
        else
            a = params(1);
            b = params(2);
        end
        
        % Stan początkowy
        if nargin < 4 || isempty(initial_state)
            x = rand();
            y = rand();
        else
            x = initial_state(1);
            y = initial_state(2);
        end
        
        % Pomiń stany przejściowe
        for i = 1:transient
            x_new = 1 - a * x^2 + b * y;
            y = x;
            x = x_new;
        end
        
        % Generowanie sekwencji
        for i = 1:length
            x_new = 1 - a * x^2 + b * y;
            y = x;
            x = x_new;
            sequence(i) = x;
        end
        
    case 'lorenz'
        % Parametry dla układu Lorenza
        if nargin < 3 || isempty(params)
            sigma = 10;
            rho = 28;
            beta = 8/3;
        else
            sigma = params(1);
            rho = params(2);
            beta = params(3);
        end
        
        % Stan początkowy
        if nargin < 4 || isempty(initial_state)
            x = rand();
            y = rand();
            z = rand();
        else
            x = initial_state(1);
            y = initial_state(2);
            z = initial_state(3);
        end
        
        % Krok całkowania (można dostosować)
        dt = 0.01;
        
        % Pomiń stany przejściowe
        for i = 1:transient
            dx = sigma * (y - x);
            dy = x * (rho - z) - y;
            dz = x * y - beta * z;
            x = x + dx * dt;
            y = y + dy * dt;
            z = z + dz * dt;
        end
        
        % Generowanie sekwencji
        for i = 1:length
            dx = sigma * (y - x);
            dy = x * (rho - z) - y;
            dz = x * y - beta * z;
            x = x + dx * dt;
            y = y + dy * dt;
            z = z + dz * dt;
            sequence(i) = x; % Można też użyć y lub z
        end
        
    case 'bernoulli'
        % Parametry dla odwzorowania Bernoulli
        if nargin < 3 || isempty(params)
            beta = 0.5; % Standardowy parametr odwzorowania Bernoulli
        else
            beta = params;
        end
        
        % Stan początkowy
        if nargin < 4 || isempty(initial_state)
            x = rand(); % Losowa wartość między 0 a 1
        else
            x = initial_state;
        end
        
        % Pomiń stany przejściowe
        for i = 1:transient
            x = mod(2 * x + beta, 1);
        end
        
        % Generowanie sekwencji
        for i = 1:length
            x = mod(2 * x + beta, 1);
            sequence(i) = x;
        end
        
    otherwise
        error('Nieznany typ sekwencji chaotycznej. Dostępne typy: logistic, tent, henon, lorenz, bernoulli');
end

end

function binary_sequence = chaoticToBinary(chaotic_sequence, threshold)
% CHAOTICTOBINARY Konwertuje sekwencję chaotyczną na sekwencję binarną
%
% Składnia:
%   binary_sequence = chaoticToBinary(chaotic_sequence, threshold)
%
% Wejście:
%   chaotic_sequence - wejściowa sekwencja chaotyczna
%   threshold - próg binaryzacji (domyślnie 0.5 lub średnia wartość sekwencji)
%
% Wyjście:
%   binary_sequence - sekwencja binarna (0 i 1)

if nargin < 2 || isempty(threshold)
    threshold = mean(chaotic_sequence);
end

binary_sequence = double(chaotic_sequence > threshold);

end

function bipolar_sequence = chaoticToBipolar(chaotic_sequence, threshold)
% CHAOTICTOBIPOLAR Konwertuje sekwencję chaotyczną na sekwencję bipolarną
%
% Składnia:
%   bipolar_sequence = chaoticToBipolar(chaotic_sequence, threshold)
%
% Wejście:
%   chaotic_sequence - wejściowa sekwencja chaotyczna
%   threshold - próg binaryzacji (domyślnie 0.5 lub średnia wartość sekwencji)
%
% Wyjście:
%   bipolar_sequence - sekwencja bipolarna (-1 i 1)

if nargin < 2 || isempty(threshold)
    threshold = mean(chaotic_sequence);
end

bipolar_sequence = 2 * double(chaotic_sequence > threshold) - 1;

end

function [Rxx, lags] = computeAutoCorrelation(sequence)
% COMPUTEAUTOCORRELATION Oblicza autokorelację sekwencji
%
% Składnia:
%   [Rxx, lags] = computeAutoCorrelation(sequence)
%
% Wejście:
%   sequence - wejściowa sekwencja
%
% Wyjście:
%   Rxx - funkcja autokorelacji
%   lags - przesunięcia (opóźnienia)

[Rxx, lags] = xcorr(sequence, 'biased');

end

function [Rxy, lags] = computeCrossCorrelation(sequence1, sequence2)
% COMPUTECROSSCORRELATION Oblicza korelację krzyżową dwóch sekwencji
%
% Składnia:
%   [Rxy, lags] = computeCrossCorrelation(sequence1, sequence2)
%
% Wejście:
%   sequence1 - pierwsza sekwencja
%   sequence2 - druga sekwencja
%
% Wyjście:
%   Rxy - funkcja korelacji krzyżowej
%   lags - przesunięcia (opóźnienia)

[Rxy, lags] = xcorr(sequence1, sequence2, 'biased');

end

function demoChaoticSequences()
% DEMOCHAOTICSEQUENCES Demonstracja generowania i analizy sekwencji chaotycznych
    
    % Parametry
    N = 1000; % Długość sekwencji
    transient = 200; % Długość przebiegu przejściowego
    
    % 1. Generowanie różnych typów sekwencji chaotycznych
    
    % Odwzorowanie logistyczne
    r_values = [3.7, 3.8, 3.9, 4.0];
    figure;
    sgtitle('Odwzorowanie Logistyczne dla różnych wartości r');
    
    for i = 1:length(r_values)
        x0 = 0.4; % Stały stan początkowy
        logistic_seq = generateChaoticSequence('logistic', N, r_values(i), x0, transient);
        
        subplot(length(r_values), 1, i);
        plot(logistic_seq(1:200));
        title(['r = ', num2str(r_values(i))]);
        xlabel('Iteracja');
        ylabel('x(n)');
        grid on;
    end
    
    % 2. Demonstracja wrażliwości na warunki początkowe (efekt motyla)
    figure;
    sgtitle('Wrażliwość na warunki początkowe (odwzorowanie logistyczne, r = 3.99)');
    
    delta = 1e-10; % Mała różnica w warunku początkowym
    x0_1 = 0.5;
    x0_2 = x0_1 + delta;
    
    logistic_seq1 = generateChaoticSequence('logistic', N, 3.99, x0_1, 0);
    logistic_seq2 = generateChaoticSequence('logistic', N, 3.99, x0_2, 0);
    
    subplot(2, 1, 1);
    plot(logistic_seq1(1:200), 'b');
    hold on;
    plot(logistic_seq2(1:200), 'r--');
    title(['Początkowa różnica = ', num2str(delta)]);
    legend('Sekwencja 1', 'Sekwencja 2');
    xlabel('Iteracja');
    ylabel('x(n)');
    grid on;
    
    subplot(2, 1, 2);
    plot(abs(logistic_seq1 - logistic_seq2));
    title('Bezwzględna różnica między sekwencjami');
    xlabel('Iteracja');
    ylabel('|x_1(n) - x_2(n)|');
    grid on;
    
    % 3. Porównanie binarnych sekwencji z różnych map chaotycznych
    figure;
    sgtitle('Binarne sekwencje z różnych map chaotycznych');
    
    % Generowanie sekwencji
    logistic_seq = generateChaoticSequence('logistic', N, 3.99, 0.5, transient);
    tent_seq = generateChaoticSequence('tent', N, 1.99, 0.5, transient);
    henon_seq = generateChaoticSequence('henon', N, [1.4, 0.3], [0.5, 0.5], transient);
    lorenz_seq = generateChaoticSequence('lorenz', N, [10, 28, 8/3], [1, 1, 1], transient);
    bernoulli_seq = generateChaoticSequence('bernoulli', N, 0.5, 0.5, transient);
    
    % Konwersja na sekwencje binarne
    bin_logistic = chaoticToBinary(logistic_seq);
    bin_tent = chaoticToBinary(tent_seq);
    bin_henon = chaoticToBinary(henon_seq);
    bin_lorenz = chaoticToBinary(lorenz_seq);
    bin_bernoulli = chaoticToBinary(bernoulli_seq);
    
    % Wyświetlanie pierwszych 100 bitów z każdej sekwencji
    subplot(5, 1, 1);
    stem(bin_logistic(1:100), 'filled');
    title('Sekwencja binarna z odwzorowania logistycznego');
    ylim([-0.1, 1.1]);
    grid on;
    
    subplot(5, 1, 2);
    stem(bin_tent(1:100), 'filled');
    title('Sekwencja binarna z odwzorowania namiotowego');
    ylim([-0.1, 1.1]);
    grid on;
    
    subplot(5, 1, 3);
    stem(bin_henon(1:100), 'filled');
    title('Sekwencja binarna z odwzorowania Hénona');
    ylim([-0.1, 1.1]);
    grid on;
    
    subplot(5, 1, 4);
    stem(bin_lorenz(1:100), 'filled');
    title('Sekwencja binarna z układu Lorenza');
    ylim([-0.1, 1.1]);
    grid on;
    
    subplot(5, 1, 5);
    stem(bin_bernoulli(1:100), 'filled');
    title('Sekwencja binarna z odwzorowania Bernoulli');
    ylim([-0.1, 1.1]);
    grid on;
    
    % 4. Analiza właściwości korelacyjnych
    figure;
    sgtitle('Właściwości korelacyjne binarnych sekwencji chaotycznych');
    
    % Konwersja na sekwencje bipolarne dla lepszych właściwości korelacyjnych
    bip_logistic = chaoticToBipolar(logistic_seq);
    bip_tent = chaoticToBipolar(tent_seq);
    
    % Obliczanie autokorelacji
    [acorr_logistic, lags_l] = computeAutoCorrelation(bip_logistic);
    [acorr_tent, lags_t] = computeAutoCorrelation(bip_tent);
    
    % Obliczanie korelacji krzyżowej
    [ccorr, lags_c] = computeCrossCorrelation(bip_logistic, bip_tent);
    
    % Wyświetlanie autokorelacji
    subplot(3, 1, 1);
    plot(lags_l, acorr_logistic);
    title('Autokorelacja sekwencji z odwzorowania logistycznego');
    xlabel('Przesunięcie');
    ylabel('Autokorelacja');
    grid on;
    
    subplot(3, 1, 2);
    plot(lags_t, acorr_tent);
    title('Autokorelacja sekwencji z odwzorowania namiotowego');
    xlabel('Przesunięcie');
    ylabel('Autokorelacja');
    grid on;
    
    % Wyświetlanie korelacji krzyżowej
    subplot(3, 1, 3);
    plot(lags_c, ccorr);
    title('Korelacja krzyżowa między sekwencjami');
    xlabel('Przesunięcie');
    ylabel('Korelacja krzyżowa');
    grid on;
    
    % 5. Atraktor układu Lorenza
    lorenz_long = generateChaoticSequence('lorenz', 10000, [10, 28, 8/3], [1, 1, 1], 1000);
    
    % Całkowanie układu Lorenza dla uzyskania wszystkich 3 współrzędnych
    sigma = 10;
    rho = 28;
    beta = 8/3;
    dt = 0.01;
    
    x = 1;
    y = 1;
    z = 1;
    
    X = zeros(1, 10000);
    Y = zeros(1, 10000);
    Z = zeros(1, 10000);
    
    % Pomiń stany przejściowe
    for i = 1:1000
        dx = sigma * (y - x);
        dy = x * (rho - z) - y;
        dz = x * y - beta * z;
        x = x + dx * dt;
        y = y + dy * dt;
        z = z + dz * dt;
    end
    
    % Generowanie sekwencji
    for i = 1:10000
        dx = sigma * (y - x);
        dy = x * (rho - z) - y;
        dz = x * y - beta * z;
        x = x + dx * dt;
        y = y + dy * dt;
        z = z + dz * dt;
        X(i) = x;
        Y(i) = y;
        Z(i) = z;
    end
    
    figure;
    plot3(X, Y, Z, 'b', 'LineWidth', 0.5);
    title('Atraktor układu Lorenza');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    grid on;
    view(27, 16);
    
    % 6. Wizualizacja odwzorowania logistycznego
    figure;
    r = 3.9;
    x0 = 0.5;
    n_iter = 200;
    x = zeros(1, n_iter+1);
    x(1) = x0;
    
    for i = 1:n_iter
        x(i+1) = r * x(i) * (1 - x(i));
    end
    
    % Wykres odwzorowania logistycznego
    subplot(2, 1, 1);
    plot(x, 'b.-');
    title(['Odwzorowanie logistyczne (r = ', num2str(r), ')']);
    xlabel('Iteracja n');
    ylabel('x(n)');
    grid on;
    
    % Wykres Cobweb (pajęczyny)
    subplot(2, 1, 2);
    t = linspace(0, 1, 1000);
    ft = r * t .* (1 - t);
    
    plot(t, ft, 'b', 'LineWidth', 1.5);  % Funkcja f(x) = r*x*(1-x)
    hold on;
    plot(t, t, 'k--');  % Linia y = x
    
    % Rysowanie diagramu pajęczyny
    for i = 1:20
        xi = x(i);
        fxi = x(i+1);
        
        % Linia pionowa od (xi, xi) do (xi, fxi)
        plot([xi, xi], [xi, fxi], 'r-');
        
        % Linia pozioma od (xi, fxi) do (fxi, fxi)
        plot([xi, fxi], [fxi, fxi], 'r-');
    end
    
    title('Diagram Cobweb (pajęczyny) dla odwzorowania logistycznego');
    xlabel('x');
    ylabel('f(x)');
    axis([0 1 0 1]);
    grid on;
end

function demoChaoticCDMA()
% DEMOCHAOTICCDMA Demonstracja zastosowania sekwencji chaotycznych w systemie CDMA
    
    % Parametry
    num_users = 4;           % Liczba użytkowników
    bits_per_user = 10;      % Liczba bitów na użytkownika
    chips_per_bit = 64;      % Liczba chipów na bit (długość sekwencji rozpraszającej)
    SNR_dB = 10;             % Stosunek sygnału do szumu w dB
    
    % Generowanie losowych bitów dla użytkowników
    data_bits = randi([0, 1], num_users, bits_per_user);
    
    % Konwersja na wartości bipolarne (-1, 1)
    data_bipolar = 2 * data_bits - 1;
    
    % Generowanie chaotycznych sekwencji rozpraszających dla każdego użytkownika
    spreading_sequences = zeros(num_users, chips_per_bit);
    
    for u = 1:num_users
        % Generujemy różne sekwencje chaotyczne dla każdego użytkownika
        % Używamy różnych wartości r i warunków początkowych dla odwzorowania logistycznego
        r = 3.9 + 0.05 * (u-1);  % Różne wartości r dla każdego użytkownika
        x0 = 0.1 + 0.2 * (u-1);  % Różne wartości początkowe
        
        chaotic_seq = generateChaoticSequence('logistic', chips_per_bit, r, x0, 100);
        spreading_sequences(u, :) = chaoticToBipolar(chaotic_seq);
    end
    
    % Rozpraszanie danych
    spread_data = zeros(num_users, bits_per_user * chips_per_bit);
    
    for u = 1:num_users
        for b = 1:bits_per_user
            chip_range = (b-1)*chips_per_bit+1 : b*chips_per_bit;
            spread_data(u, chip_range) = data_bipolar(u, b) * spreading_sequences(u, :);
        end
    end
    
    % Sumowanie sygnałów wszystkich użytkowników (kanał CDMA)
    composite_signal = sum(spread_data, 1);
    
    % Dodanie szumu gaussowskiego
    SNR_linear = 10^(SNR_dB/10);
    signal_power = var(composite_signal);
    noise_power = signal_power / SNR_linear;
    noise = sqrt(noise_power) * randn(size(composite_signal));
    
    received_signal = composite_signal + noise;
    
    % Odbiornik - derozpraszanie dla każdego użytkownika
    recovered_data = zeros(num_users, bits_per_user);
    
    for u = 1:num_users
        for b = 1:bits_per_user
            chip_range = (b-1)*chips_per_bit+1 : b*chips_per_bit;
            corr_value = received_signal(chip_range) * spreading_sequences(u, :)' / chips_per_bit;
            recovered_data(u, b) = sign(corr_value);  % Decyzja twarda
        end
    end
    
    % Obliczenie BER (Bit Error Rate) dla każdego użytkownika
    errors = sum(sum(data_bipolar ~= recovered_data)) / 2;  % Dzielimy przez 2, bo -1 -> 1 i 1 -> -1 to dwa błędy
    total_bits = num_users * bits_per_user;
    BER = errors / total_bits;
    
    % Wyświetlanie wyników
    fprintf('Symulacja systemu CDMA z sekwencjami chaotycznymi:\n');
    fprintf('Liczba użytkowników: %d\n', num_users);
    fprintf('Liczba bitów na użytkownika: %d\n', bits_per_user);
    fprintf('Długość sekwencji rozpraszającej (chips/bit): %d\n', chips_per_bit);
    fprintf('SNR: %g dB\n', SNR_dB);
    fprintf('Całkowity BER: %g\n', BER);
    
    % Wizualizacja
    figure;
    sgtitle('Symulacja systemu CDMA z sekwencjami chaotycznymi');
    
    % Wyświetl sekwencje rozpraszające
    subplot(4, 1, 1);
    imagesc(spreading_sequences);
    colormap(gray);
    title('Sekwencje rozpraszające dla użytkowników');
    xlabel('Indeks chipu');
    ylabel('Użytkownik');
    colorbar;
    
    % Wyświetl dane przed rozproszeniem
    subplot(4, 1, 2);
    imagesc(data_bipolar);
    title('Dane użytkowników (dane bipolarne)');
    xlabel('Indeks bitu');
    ylabel('Użytkownik');
    colorbar;
    
    % Wyświetl sygnał zbiorczy (z kanału)
    subplot(4, 1, 3);
    plot(received_signal(1:min(500, length(received_signal))));
    title('Fragment odebranego sygnału zbiorczego');
    xlabel('Indeks próbki');
    ylabel('Amplituda');
    grid on;
    
    % Wyświetl odzyskane dane
    subplot(4, 1, 4);
    imagesc(recovered_data);
    title('Odzyskane dane użytkowników');
    xlabel('Indeks bitu');
    ylabel('Użytkownik');
    colorbar;
    
    % Korelacja między sekwencjami rozpraszającymi
    figure;
    sgtitle('Właściwości korelacyjne sekwencji rozpraszających');
    
    % Autokorelacja pierwszej sekwencji
    [acorr, lags_a] = computeAutoCorrelation(spreading_sequences(1, :));
    
    subplot(2, 1, 1);
    plot(lags_a, acorr);
    title('Autokorelacja pierwszej sekwencji rozpraszającej');
    xlabel('Przesunięcie');
    ylabel('Autokorelacja');
    grid on;
    
    % Korelacja krzyżowa między pierwszą a drugą sekwencją
    [ccorr, lags_c] = computeCrossCorrelation(spreading_sequences(1, :), spreading_sequences(2, :));
    
    subplot(2, 1, 2);
    plot(lags_c, ccorr);
    title('Korelacja krzyżowa: sekwencja 1 i 2');
    xlabel('Przesunięcie');
    ylabel('Korelacja krzyżowa');
    grid on;
end