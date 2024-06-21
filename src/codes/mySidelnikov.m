% 
% % Generowanie 5 sekwencji Sidel'nikova o długości 100 bitów dla p=101
% [sequences, primitive_elements] = mySidelnikov(101, 1, 5, 100);
% 
% % Wyświetlenie pierwszej sekwencji
% figure;
% stem(sequences(1,:));
% title('Sekwencja Sidel''nikova');
% xlabel('Indeks');
% ylabel('Wartość bitu');


function sidelnikovSequences = mySidelnikov(q, d, numSeq)
% sidelnikovGenerator generuje rodzinę sekwencji Sidelnikova.
%
%   sidelnikovSequences = sidelnikovGenerator(q, d, numSeq)
%
%   Wejście:
%       q      - rozmiar ciała GF(q) (q musi być potęgą liczby pierwszej)
%       d      - dzielnik liczby (q-1); element beta będzie miał rząd równy d
%       numSeq - liczba sekwencji do wygenerowania (maksymalnie q różnych przesunięć gamma)
%
%   Wyjście:
%       sidelnikovSequences - macierz o wymiarach numSeq x (q-1), gdzie każdy
%                             wiersz to jedna sekwencja Sidelnikova
%
% Przykład użycia:
%   q = 7;   % Rozmiar ciała GF(7)
%   d = 3;   % Dzielnik (q-1), czyli 6
%   numSeq = 4;  % Liczba sekwencji do wygenerowania
%   sequences = sidelnikovGenerator(q, d, numSeq);
%   disp(sequences);
%
% Wymagane: Communications Toolbox (funkcja gf, primpoly) dla GF(2^m)

    % Sprawdź, czy q jest potęgą liczby pierwszej
    if ~isPrimePower(q)
        error('q musi być potęgą liczby pierwszej');
    end
    
    % Określenie, czy pracujemy w GF(p) czy GF(p^m)
    if isprime(q)
        p = q;
        m = 1;
    else
        % Znajdź p i m, takie że q = p^m
        factors = factor(q);
        p = factors(1);
        if ~all(factors == p)
            error('q musi być potęgą liczby pierwszej');
        end
        m = length(factors);
    end
    
    L = q - 1; % długość sekwencji
    
    % Sprawdź, czy d jest dzielnikiem (q-1)
    if mod(L, d) ~= 0
        error('d musi być dzielnikiem (q-1)');
    end
    
    % Przygotujemy tablicę potęg generatora
    powers = zeros(L, 1);
    log_table = zeros(q, 1);
    
    % Dla GF(p) gdzie p jest liczbą pierwszą, używamy generatora
    if m == 1
        % Znajdź generator (element pierwotny) dla GF(p)
        g = findPrimitiveElementForPrime(p);
        % Oblicz wszystkie potęgi g^i mod p
        power_val = 1;
        for i = 1:L
            power_val = mod(power_val * g, p);
            powers(i) = power_val;
            log_table(power_val + 1) = i; % +1 bo indeksowanie od 1
        end
        
        % Ustal beta jako element o rządzie d: beta = g^((p-1)/d)
        beta_power = L / d;
        beta = mod(g^beta_power, p);
        
        % Przygotuj macierz wynikową
        sidelnikovSequences = zeros(numSeq, L);
        
        % Sprawdź, czy numSeq nie przekracza liczby elementów w GF(q)
        if numSeq > q
            error('Liczba żądanych sekwencji przekracza liczbę elementów GF(%d).', q);
        end
        
        % Dla każdej wartości gamma generujemy sekwencję
        for j = 1:numSeq
            gamma = j - 1; % gamma przyjmuje wartości 0, 1, 2, ..., numSeq-1
            seq = zeros(1, L);
            
            for i = 1:L
                % a = g^i + gamma (operacje w GF(p))
                a = mod(powers(i) + gamma, p);
                
                if a == 0
                    seq(i) = 0;
                else
                    % logarytm: find k such that g^k = a
                    k = log_table(a + 1);
                    
                    % Oblicz logarytm względem beta (rzędu d)
                    r = mod(round(k * d / L), d);
                    seq(i) = r;
                end
            end
            
            sidelnikovSequences(j, :) = seq;
        end
    else
        % Dla GF(2^m) z m > 1, używamy Communications Toolbox
        if p > 2
            error('Obsługa GF(p^m) dla p > 2 nie jest obecnie zaimplementowana.');
        end
        
        if ~license('test', 'Communication_Toolbox')
            error('Do obsługi GF(2^m) wymagany jest Communications Toolbox.');
        end
        
        prim_poly = primpoly(m);
        
        % Tworzymy elementy ciała GF(2^m)
        field_elements = gf([0:2^m-1]', m, prim_poly);
        
        % Znajdź element pierwotny (generator) ciała GF(2^m)
        generator_idx = findPrimitiveElementIdx(field_elements, m, prim_poly);
        alpha = field_elements(generator_idx + 1);
        
        % Przygotuj tablicę potęg i logarytmów
        powers = zeros(L, 1);
        log_table = zeros(2^m, 1);
        
        % Oblicz alpha^i dla i=1,2,...,L
        power = alpha;
        for i = 1:L
            powers(i) = gf2idx(power, m);
            log_table(powers(i) + 1) = i;
            power = power * alpha;
        end
        
        % beta = alpha^((2^m-1)/d) - element rzędu d
        beta_power = round(L / d);
        beta = alpha^beta_power;
        
        % Przygotuj macierz wynikową
        sidelnikovSequences = zeros(numSeq, L);
        
        % Dla każdej wartości gamma (przesunięcia)
        for j = 1:numSeq
            gamma_idx = j - 1;
            gamma = field_elements(gamma_idx + 1);
            
            seq = zeros(1, L);
            for i = 1:L
                % a = alpha^i + gamma
                ai = alpha^i;
                a = ai + gamma;
                
                % Sprawdź, czy a = 0
                if a == 0
                    seq(i) = 0;
                else
                    % Znajdź k takie, że a = alpha^k
                    a_idx = gf2idx(a, m);
                    k = log_table(a_idx + 1);
                    
                    % Oblicz logarytm względem beta (element rzędu d)
                    r = mod(round(k * d / L), d);
                    seq(i) = r;
                end
            end
            
            sidelnikovSequences(j, :) = seq;
        end
    end
end

% Funkcja sprawdzająca, czy liczba jest potęgą liczby pierwszej
function result = isPrimePower(n)
    if n <= 1
        result = false;
        return;
    end
    
    factors = factor(n);
    first_factor = factors(1);
    result = all(factors == first_factor);
end

% Funkcja znajdująca element pierwotny (generator) dla GF(p)
function g = findPrimitiveElementForPrime(p)
    % Dla małych pól można po prostu przeiterować przez możliwe generatory
    for candidate = 2:p-1
        % Sprawdź, czy candidate jest generatorem
        powers = zeros(1, p-1);
        power_val = 1;
        is_generator = true;
        
        for i = 1:p-2
            power_val = mod(power_val * candidate, p);
            powers(i) = power_val;
            
            % Jeśli przed dotarciem do p-1 pojawi się 1, to nie jest generator
            if power_val == 1
                is_generator = false;
                break;
            end
        end
        
        % Ostatnia potęga powinna być 1 dla generatora
        power_val = mod(power_val * candidate, p);
        if power_val == 1 && is_generator
            g = candidate;
            return;
        end
    end
    
    error('Nie znaleziono elementu pierwotnego dla GF(%d).', p);
end

% Funkcja do znalezienia indeksu elementu pierwotnego w GF(2^m)
function idx = findPrimitiveElementIdx(field_elements, m, prim_poly)
    q = 2^m;
    one = field_elements(2); % indeks 1, ponieważ indeksujemy od 0
    
    % Sprawdź wszystkie niezerowe elementy
    for idx = 1:q-1
        candidate = field_elements(idx+1);
        is_generator = true;
        
        power = candidate;
        for ord = 1:q-2
            if power == one
                is_generator = false;
                break;
            end
            power = power * candidate;
        end
        
        % Ostatnia potęga powinna być 1
        if is_generator && power == one
            return;
        end
    end
    
    error('Nie znaleziono elementu pierwotnego w GF(2^%d).', m);
end

% Funkcja konwersji elementu GF(2^m) na indeks (0...2^m-1)
function idx = gf2idx(element, m)
    % Konwertuj wektor binarny na liczbę
    bin_vec = element.x;
    idx = 0;
    for i = 1:length(bin_vec)
        idx = idx + bin_vec(i) * 2^(i-1);
    end
end