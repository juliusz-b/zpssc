%% Zadanie_DFE_test.m
% Test DFE/SIC do kompensacji sidelobes i ghostów w systemie CDM-FBG
% Podejście: equalizacja sygnału kanałowego PRZED korelacją
%
% Idea: kanał FBG = filtr z wieloma odbiciami (jak ISI w telekomunikacji)
%       DFE uczy się na znanej sekwencji treningowej (kod Kasami)
%       i kompensuje wielokrotne odbicia
%
% Projekt: PN/01/0321/2022
% Data: 2025-12-15

clear all; close all; clc;
cd('C:\Users\jbojarczuk\OneDrive - Politechnika Warszawska\PRACA\Insytut Telekomunikacji\Perly Nauki 2022\zpssc');
AddAllSubfolders();

SILENT_MODE = true;
neff = 1.447;

%% ===== KROK 1: SYMULACJA BAZOWA =====
% Prosty scenariusz: 5 siatek, p=8, domyślne parametry
fprintf('=== KROK 1: Symulacja bazowa ===\n');

scenariusz = 'kody2'; kod = 'kasami'; mode = 'unipolar';
p = 8; Nb = 8; Fsample = 4e9; K = 1;

laser.lasers_mode = 'number'; laser.step = 16;
laser.wavelength_range = [1534.5, 1536.5];
laser.FWHM = 2.4; laser.power = 10^(-12/10)*1e-3; laser.shape = 'lorentz';

pd.A = 0.8; pd.BW = 250e6; pd.gain = 1; pd.RL = 50; pd.Idark = 0.05e-12;
fiber.alpha = 1;

fbg.N_s = 5;
fbg.D = 400:20:480;
fbg.periods = 530.73e-9 * ones(1, fbg.N_s);
fbg.deltaneffs = 0.3e-4 * ones(1, fbg.N_s);
fbg.grating_lengths = 3e-3 * ones(1, fbg.N_s);
fbg.model = 'tanh';
fbg.max_bounces = 1;
fbg.lambda_shifts = zeros(1, fbg.N_s);

noise.type = 'nep'; noise.SNR = 10; noise.NEP = 15e-12; noise.Fn = 1;
initial_filter.type = 'none'; initial_filter.alpha = 0;
modulation.method = 'amssb'; modulation.A0 = 0.1;
modulation.Fn_s = (20:90:(20+90*(laser.step-1)))*1e6;

% --- Pipeline ręcznie (skrypt, nie funkcja!) ---
L = 2^p - 1;
lambdas = generateLambdas(laser);
[simulation_range, lasers_spectra] = generateLasers(laser);
[R_s, R_simulated] = generateAllGratings(lambdas, fbg, simulation_range);
U = length(lambdas);
[D_s, L_max, time_domain, length_domain] = generateValues(L, Fsample, Nb, fbg.D, U, fbg.N_s, false, 1);
[BipUni, power] = selectMode(mode);
kody = BipUni(genSpreadCodes(L, U, kod));
dane_upsamplowane = myUpsample(kody, Nb);
[filter_data] = generateInitialFilter(initial_filter, Fsample, Nb, 16);
dane_wstepnie_filtrowane = dane_upsamplowane;

% Odbicia od siatek
dane_z_siatkami = addAllGratings(dane_wstepnie_filtrowane, fbg.N_s, D_s, R_s, L_max, ...
    R_simulated, lasers_spectra, simulation_range, fbg, fiber);

% Kanał kody2: sumowanie + fotodetektor + szum + lowpass
dane_w_kanale = sum(dane_z_siatkami, 1) * pd.A;
dane_w_kanale = noiseAndMeanDataInChannel(dane_w_kanale, K, laser, noise, pd).^power;
b_lpf = fir1(64, pd.BW / (Fsample/2));
for i = 1:size(dane_w_kanale,1)
    dane_w_kanale(i,:) = filter(b_lpf, 1, dane_w_kanale(i,:));
end
dane_w_kanale = dane_w_kanale * pd.gain;

fprintf('Sygnał kanałowy: [%d x %d]\n', size(dane_w_kanale));
fprintf('Pozycje siatek (samples): %s\n', mat2str(D_s(1:fbg.N_s)));

%% ===== KROK 2: KORELACJA BEZ EQUALIZACJI (BASELINE) =====
fprintf('\n=== KROK 2: Korelacja baseline ===\n');

xcv_baseline = calcXcov(dane_w_kanale, dane_upsamplowane, true);
xcv_baseline_bez = xcv_baseline - movmean(xcv_baseline', D_s(2)-D_s(1))';

% Detekcja i MAE
[~, R_rec_base] = findGratings(xcv_baseline_bez, dane_upsamplowane, length_domain, fbg.N_s, D_s);
R_rec_base = R_rec_base.^(1/power) / laser.power / pd.A;

lB_org = zeros(1, fbg.N_s);
lB_rec_base = zeros(1, fbg.N_s);
for i = 1:fbg.N_s
    lB_org(i) = 2*neff*fbg.periods(i)*(1+fbg.deltaneffs(i)/neff)*1e9;
    y = R_rec_base(i,:); mx = max(y)*0.3; I = y>=mx;
    if any(I); lB_rec_base(i) = sum(y(I).*lambdas(I)*1e9)/sum(y(I)); else; lB_rec_base(i)=NaN; end
end
MAE_baseline = mean(abs(lB_rec_base - lB_org)) * 1000;
fprintf('MAE baseline: %.1f pm\n', MAE_baseline);

%% ===== KROK 3: SIC (Successive Interference Cancellation) =====
% Podejście: po korelacji, znajdź najsilniejszy peak, zregeneruj jego
% wkład do sygnału korelacyjnego, odejmij, powtórz
fprintf('\n=== KROK 3: SIC na korelogramie ===\n');

% Dla każdego kanału lambda osobno
xcv_sic = xcv_baseline_bez; % kopia

for iter = 1:fbg.N_s
    % Znajdź najsilniejszy peak
    [peak_val, peak_idx] = max(max(abs(xcv_sic), [], 1));

    % Zregeneruj wkład tego peaka: autokorelacja kodu × amplituda
    % Autokorelacja kodu referencyjnego (dla kanału 1)
    ref_autocorr = xcorr(dane_upsamplowane(1,:), dane_upsamplowane(1,:));
    ref_autocorr = ref_autocorr(end/2:end); % positive lags
    ref_autocorr = ref_autocorr / max(ref_autocorr); % normalizuj

    % Dla każdego kanału lambda: odejmij wkład tego peaka
    for ch = 1:size(xcv_sic, 1)
        amp = xcv_sic(ch, peak_idx); % amplituda na pozycji peaka

        % Regeneracja: przesuń autokorelację na pozycję peaka
        regen = zeros(1, size(xcv_sic, 2));
        n_copy = min(length(ref_autocorr), size(xcv_sic,2) - peak_idx + 1);
        regen(peak_idx:peak_idx+n_copy-1) = amp * ref_autocorr(1:n_copy);
        % Lustrzane odbicie (lewa strona peaka)
        n_left = min(length(ref_autocorr)-1, peak_idx-1);
        regen(peak_idx-n_left:peak_idx-1) = amp * ref_autocorr(n_left+1:-1:2);

        xcv_sic(ch,:) = xcv_sic(ch,:) - regen;
    end

    fprintf('  Iteracja %d: usunięto peak na pozycji %.1f m (amp=%.2e)\n', ...
        iter, length_domain(peak_idx), peak_val);
end

% Po SIC: dodaj z powrotem peaki (SIC usunął wszystkie - chcemy "czyste" peaki)
% Alternatywa: użyj oryginalnych pozycji i amplitud z SIC-cleaned signal
% Na razie: po prostu porównaj korelogram przed/po SIC

%% ===== KROK 4: DFE na sygnale kanałowym =====
fprintf('\n=== KROK 4: DFE na sygnale kanałowym ===\n');

% Podejście: DFE equalizuje sygnał kanałowy PRZED korelacją
% Sekwencja treningowa = znany nadany kod (kanał 1)
% Sygnał odbierany = dane_w_kanale (po fotodetektorze)
%
% DFE uczy się "kanału" (siatki FBG + multi-bounce) i kompensuje go

% Weźmy jeden kanał (sumowany - kody2)
rx_signal = dane_w_kanale(1, :); % sygnał odebrany
tx_signal = dane_upsamplowane(1, :); % znany nadany kod (kanał 1)

% Trzeba dopasować długości
n_sig = min(length(rx_signal), length(tx_signal));

% DFE parametry
nForwardTaps = 7;   % tapy forward
nFeedbackTaps = 5;  % tapy feedback

try
    % Użyj comm.DecisionFeedbackEqualizer
    dfe = comm.DecisionFeedbackEqualizer( ...
        'Algorithm', 'LMS', ...
        'NumForwardTaps', nForwardTaps, ...
        'NumFeedbackTaps', nFeedbackTaps, ...
        'StepSize', 0.001, ...
        'ReferenceTap', ceil(nForwardTaps/2));

    % Training: znamy tx_signal
    [eq_signal, err_signal] = dfe(rx_signal(1:n_sig).', tx_signal(1:n_sig).');
    eq_signal = eq_signal.';

    fprintf('DFE wyequalizowany sygnał: [1 x %d]\n', length(eq_signal));
    fprintf('Średni błąd DFE: %.4e\n', mean(abs(err_signal)));

    % Korelacja po DFE
    % Musimy skorelować z ORYGINALNYM kodem (nie zequalizowanym)
    % Dopasuj długości
    eq_padded = [eq_signal, zeros(1, size(dane_w_kanale,2) - length(eq_signal))];

    xcv_dfe = calcXcov(eq_padded, dane_upsamplowane, true);
    xcv_dfe_bez = xcv_dfe - movmean(xcv_dfe', D_s(2)-D_s(1))';

    has_dfe = true;
    fprintf('Korelacja po DFE gotowa.\n');

catch ME
    fprintf('DFE error: %s\n', ME.message);
    fprintf('Próbuję alternatywne podejście...\n');
    has_dfe = false;
end

%% ===== KROK 5: Alternatywa - prosty equalizer ręczny =====
% Jeśli comm.DecisionFeedbackEqualizer nie działa lub daje złe wyniki,
% implementujemy prostą wersję SIC na korelogramie
fprintf('\n=== KROK 5: SIC ręczny na korelogramie ===\n');

% Klasyczny SIC: iteracyjne odejmowanie
xcv_clean = xcv_baseline; % PRZED usunięciem obwiedni
detected_positions = zeros(1, fbg.N_s);
detected_amplitudes = zeros(size(xcv_baseline, 1), fbg.N_s);

for iter = 1:fbg.N_s
    % Usuń obwiednię
    xcv_temp = xcv_clean - movmean(xcv_clean', D_s(2)-D_s(1))';

    % Znajdź najsilniejszy peak (po wszystkich kanałach)
    [~, peak_idx] = max(sum(abs(xcv_temp), 1));
    detected_positions(iter) = peak_idx;

    % Zapamiętaj amplitudę per kanał
    for ch = 1:size(xcv_clean, 1)
        detected_amplitudes(ch, iter) = xcv_temp(ch, peak_idx);
    end

    % Zregeneruj i odejmij wkład tego peaka z ORYGINALNEGO korelogramu
    for ch = 1:size(xcv_clean, 1)
        % Autokorelacja kodu dla tego kanału
        ac = xcorr(dane_upsamplowane(ch,:));
        ac = ac(ceil(length(ac)/2):end);
        ac_norm = ac / ac(1); % normalizuj: peak = 1

        amp = detected_amplitudes(ch, iter);

        % Zregeneruj sygnał korelacyjny od tej siatki
        regen = zeros(1, size(xcv_clean, 2));
        n_right = min(length(ac_norm), size(xcv_clean,2) - peak_idx);
        regen(peak_idx:peak_idx+n_right-1) = amp * ac_norm(1:n_right);
        n_left = min(length(ac_norm)-1, peak_idx-1);
        regen(peak_idx-n_left:peak_idx-1) = amp * ac_norm(n_left+1:-1:2);

        xcv_clean(ch,:) = xcv_clean(ch,:) - regen;
    end

    fprintf('  SIC iter %d: peak na %.1f m\n', iter, length_domain(peak_idx));
end

% Po SIC: zrekonstruuj czysty korelogram z wykrytymi peakami
% (oryginalne pozycje + amplitudy, ale bez sidelobes)
xcv_reconstructed = zeros(size(xcv_baseline));
for iter = 1:fbg.N_s
    pidx = detected_positions(iter);
    for ch = 1:size(xcv_reconstructed, 1)
        xcv_reconstructed(ch, pidx) = detected_amplitudes(ch, iter);
    end
end

% Alternatywa: użyj "oczyszczonego" korelogramu (residuum po SIC + peaki)
xcv_sic_final = xcv_clean; % residuum
% Dodaj peaki z powrotem (czyste, bez sidelobes)
for iter = 1:fbg.N_s
    pidx = detected_positions(iter);
    for ch = 1:size(xcv_sic_final, 1)
        xcv_sic_final(ch, pidx) = xcv_sic_final(ch, pidx) + detected_amplitudes(ch, iter);
    end
end
xcv_sic_final_bez = xcv_sic_final - movmean(xcv_sic_final', D_s(2)-D_s(1))';

% Detekcja MAE po SIC
[~, R_rec_sic] = findGratings(xcv_sic_final_bez, dane_upsamplowane, length_domain, fbg.N_s, D_s);
R_rec_sic = R_rec_sic.^(1/power) / laser.power / pd.A;

lB_rec_sic = zeros(1, fbg.N_s);
for i = 1:fbg.N_s
    y = R_rec_sic(i,:); mx = max(y)*0.3; I = y>=mx;
    if any(I); lB_rec_sic(i) = sum(y(I).*lambdas(I)*1e9)/sum(y(I)); else; lB_rec_sic(i)=NaN; end
end
MAE_sic = mean(abs(lB_rec_sic - lB_org)) * 1000;
fprintf('\nMAE po SIC: %.1f pm (vs baseline %.1f pm)\n', MAE_sic, MAE_baseline);
fprintf('Poprawa SIC: %.1f%%\n', (MAE_baseline - MAE_sic)/MAE_baseline*100);

%% ===== KROK 6: WYKRESY =====
fprintf('\n=== Generowanie wykresów ===\n');

ch = 1;
figure('Position', [50 50 1400 800], 'theme', 'light', 'Name', 'SIC/DFE Test');

% Baseline vs SIC
subplot(3,1,1);
plot(length_domain, xcv_baseline_bez(ch,:), 'b-', 'LineWidth', 1.5); hold on;
plot(length_domain, xcv_sic_final_bez(ch,:), 'r-', 'LineWidth', 1.2);
for j = 1:fbg.N_s
    xline(length_domain(D_s(j)), '--k', sprintf('FBG%d', j), 'LineWidth', 1, 'FontSize', 9);
end
xlim([min(fbg.D)-30, max(fbg.D)+30]);
xlabel('Pozycja [m]'); ylabel('Amplituda');
title(sprintf('Korelacja: Baseline (MAE=%.1f pm) vs SIC (MAE=%.1f pm)', MAE_baseline, MAE_sic));
legend('Baseline', 'Po SIC');
grid on; set(gca, 'FontSize', 12);

% Residuum SIC (co zostało po odejmowaniu)
subplot(3,1,2);
xcv_residuum = xcv_clean - movmean(xcv_clean', D_s(2)-D_s(1))';
plot(length_domain, xcv_residuum(ch,:), 'k-', 'LineWidth', 1);
xlim([min(fbg.D)-30, max(fbg.D)+30]);
xlabel('Pozycja [m]'); ylabel('Amplituda');
title('Residuum po SIC (sidelobes + szum)');
grid on; set(gca, 'FontSize', 12);

% Porównanie reflektancji
subplot(3,1,3);
lmb_nm = lambdas * 1e9;
plot(lmb_nm, R_rec_base(1,:), 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on;
plot(lmb_nm, R_rec_sic(1,:), 'rs-', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
xlabel('\lambda [nm]'); ylabel('Reflektancja');
title('Odtworzone widmo FBG1: Baseline vs SIC');
legend('Baseline', 'Po SIC');
grid on; set(gca, 'FontSize', 12);

% DFE jeśli dostępny
if has_dfe
    figure('Position', [50 50 1200 400], 'theme', 'light', 'Name', 'DFE result');
    plot(length_domain, xcv_baseline_bez(ch,:), 'b-', 'LineWidth', 1.5); hold on;
    plot(length_domain(1:size(xcv_dfe_bez,2)), xcv_dfe_bez(ch,:), 'r-', 'LineWidth', 1.2);
    xlim([min(fbg.D)-30, max(fbg.D)+30]);
    xlabel('Pozycja [m]'); ylabel('Amplituda');
    title('Korelacja: Baseline vs po DFE');
    legend('Baseline', 'Po DFE');
    grid on; set(gca, 'FontSize', 12);
end

%% ===== PODSUMOWANIE =====
fprintf('\n========================================\n');
fprintf('  PODSUMOWANIE SIC/DFE TEST\n');
fprintf('========================================\n');
fprintf('  MAE baseline: %.1f pm\n', MAE_baseline);
fprintf('  MAE po SIC:   %.1f pm\n', MAE_sic);
fprintf('  Poprawa SIC:  %.1f%%\n', (MAE_baseline - MAE_sic)/MAE_baseline*100);
if has_dfe
    % Oblicz MAE po DFE
    [~, R_rec_dfe] = findGratings(xcv_dfe_bez, dane_upsamplowane, length_domain, fbg.N_s, D_s);
    R_rec_dfe = R_rec_dfe.^(1/power) / laser.power / pd.A;
    lB_rec_dfe = zeros(1, fbg.N_s);
    for i = 1:fbg.N_s
        y = R_rec_dfe(i,:); mx = max(y)*0.3; I = y>=mx;
        if any(I); lB_rec_dfe(i) = sum(y(I).*lambdas(I)*1e9)/sum(y(I)); else; lB_rec_dfe(i)=NaN; end
    end
    MAE_dfe = mean(abs(lB_rec_dfe - lB_org)) * 1000;
    fprintf('  MAE po DFE:   %.1f pm\n', MAE_dfe);
    fprintf('  Poprawa DFE:  %.1f%%\n', (MAE_baseline - MAE_dfe)/MAE_baseline*100);
end
fprintf('========================================\n');
