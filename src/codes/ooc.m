function codes = ooc(n,w,lambda)


%% Generowanie wszystkich kandydatów - wektory binarne o wadze w
combIndices = nchoosek(1:n, w);    % każda kombinacja pozycji jedynek
numCandidates = size(combIndices, 1);
codes = [];   % zbiór zaakceptowanych kodów (każdy wiersz to kod)


%% Główna pętla generująca kody OOC (algorytm zachłanny)
for i = 1:numCandidates
    % Tworzymy kandydat - wektor binarny długości n z jedynkami na pozycjach podanych przez combIndices(i,:)
    candidate = zeros(1, n);
    candidate(combIndices(i,:)) = 1;
    
    % Sprawdzamy autokorelację dla tego kandydata
    if ~checkAutocorrelation(candidate, lambda)
        continue;  % kandydat nie spełnia warunku autokorelacji
    end
    
    % Sprawdzamy krzyżową korelację ze wszystkimi już zaakceptowanymi kodami
    acceptCandidate = true;
    for j = 1:size(codes, 1)
        if ~checkCrossCorrelation(candidate, codes(j, :), lambda)
            acceptCandidate = false;
            break;
        end
    end
    
    % Jeśli kandydat przeszedł oba testy, dodajemy go do zbioru
    if acceptCandidate
        codes = [codes; candidate];
    end
end





end


%% Funkcje pomocnicze jako funkcje lokalne (można je umieścić w oddzielnych plikach)
% Sprawdza, czy wektor x spełnia warunek autokorelacji
function ok = checkAutocorrelation(x, lambda)
    n = length(x);
    ok = true;
    for tau = 1:n-1   % pomijamy tau = 0, bo wtedy suma = w
        ac = dot(x, circshift(x, [0, tau]));
        if ac > lambda
            ok = false;
            return;
        end
    end
end

% Sprawdza, czy wektory x i y spełniają warunek krzyżowej korelacji
function ok = checkCrossCorrelation(x, y, lambda)
    n = length(x);
    ok = true;
    for tau = 0:n-1
        cc = dot(x, circshift(y, [0, tau]));
        if cc > lambda
            ok = false;
            return;
        end
    end
end