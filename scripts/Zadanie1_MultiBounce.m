%% Zadanie1_MultiBounce.m
% Demonstracja modelu wielokrotnych odbić (multi-bounce / ghosting)
% Porównanie max_bounces=1 vs max_bounces=3
%
% Wyniki:
%   - Ghosty nakładające się na primary peaks (FBG3: +6.8%, FBG4: +15%, FBG5: +21%)
%   - Ghosty poza siatkami (700, 800, 900, 1000 m)
%   - Im dalsze siatki, tym więcej ghostów się na nie nakłada
%
% Projekt: PN/01/0321/2022

clear all; close all; clc;
cd('C:\Users\jbojarczuk\OneDrive - Politechnika Warszawska\PRACA\Insytut Telekomunikacji\Perly Nauki 2022\zpssc');
AddAllSubfolders();

%% ===== PARAMETRY =====
params = defaultParams();
params.show_plots = false;

% Silne siatki (R≈30%), spacing 100m - najlepiej widoczne ghosty
params.fbg.N_s = 5;
params.fbg.D = 200:100:600;
params.fbg.periods = 530.7300e-9 * ones(1, 5);
params.fbg.deltaneffs = 1e-4 * ones(1, 5);  % R_peak ≈ 30%
params.fbg.grating_lengths = 3e-3 * ones(1, 5);
params.fbg.lambda_shifts = zeros(1, 5);

%% ===== SYMULACJE =====
% max_bounces=1
params.fbg.max_bounces = 1;
fprintf('Obliczanie max_bounces=1... '); tic;
out1 = runSimulation(params);
fprintf('MAE=%.1f pm, t=%.1fs\n', out1.MAE, toc);

% max_bounces=3
params.fbg.max_bounces = 3;
fprintf('Obliczanie max_bounces=3... '); tic;
out3 = runSimulation(params);
fprintf('MAE=%.1f pm, t=%.1fs\n', out3.MAE, toc);

%% ===== ANALIZA =====
ch = 1;
peak_window = 30;
D_s = out3.D_s;
length_domain = out3.length_domain;
dp1 = out1.dane_bez_obwiedni;
dp3 = out3.dane_bez_obwiedni;

% Amplitudy primary peaks
amp_1b = zeros(1, params.fbg.N_s);
amp_3b = zeros(1, params.fbg.N_s);
for j = 1:params.fbg.N_s
    win = max(1, D_s(j)-peak_window):min(size(dp1,2), D_s(j)+peak_window);
    amp_1b(j) = max(dp1(ch, win));
    amp_3b(j) = max(dp3(ch, win));
end
change_pct = (amp_3b - amp_1b) ./ amp_1b * 100;

% Ghosty per siatka
ghost_counts = zeros(1, params.fbg.N_s);
ghost_pos_all = [];
for j = 1:params.fbg.N_s
    for k = 1:j-1
        for m = k+1:params.fbg.N_s
            gp = params.fbg.D(j) - params.fbg.D(k) + params.fbg.D(m);
            ghost_pos_all(end+1) = gp;
            for n = 1:params.fbg.N_s
                if abs(gp - params.fbg.D(n)) < 1
                    ghost_counts(n) = ghost_counts(n) + 1;
                end
            end
        end
    end
end
ghost_beyond = setdiff(unique(round(ghost_pos_all)), round(params.fbg.D));

%% ===== WYNIKI TEKSTOWE =====
fprintf('\n========================================\n');
fprintf('  WYNIKI MULTI-BOUNCE (R≈30%%, spacing=100m)\n');
fprintf('========================================\n');
fprintf('%-8s %-10s %-12s %-12s %-10s\n', 'Siatka', 'Ghosty', 'Amp(1b)', 'Amp(3b)', 'Zmiana%');
for j = 1:params.fbg.N_s
    fprintf('FBG%-5d %-10d %-12.4e %-12.4e %-+10.2f\n', j, ghost_counts(j), amp_1b(j), amp_3b(j), change_pct(j));
end
fprintf('\nGhosty poza siatkami: %s m\n', mat2str(ghost_beyond));

%% ===== WYKRESY =====
colors_g = [0 0.6 0; 0 0.6 0; 1 0.5 0; 1 0.2 0; 0.8 0 0];

figure('Position', [50 50 1500 900], 'theme', 'light', 'Name', 'Zad.1: Multi-bounce porównanie');

subplot(3,1,1);
plot(length_domain, dp1(ch,:), 'b-', 'LineWidth', 1.5); hold on;
plot(length_domain, dp3(ch,:), 'r-', 'LineWidth', 1.2);
for j = 1:params.fbg.N_s
    xline(length_domain(D_s(j)), '-', sprintf('FBG%d', j), ...
        'LineWidth', 1.2, 'Color', [0.3 0.3 0.3], 'FontSize', 10, 'LabelOrientation', 'horizontal');
end
for g = 1:length(ghost_beyond)
    xline(ghost_beyond(g), ':', sprintf('G%.0fm', ghost_beyond(g)), ...
        'LineWidth', 1.2, 'Color', [0.8 0 0.6], 'FontSize', 9, 'LabelOrientation', 'horizontal');
end
xlabel('Pozycja [m]'); ylabel('Amplituda');
title(sprintf('Korelacja krzyżowa \\lambda_1 | \\Deltan_{eff}=1e-4, R\\approx30%%, spacing=100m'));
xlim([100, max(ghost_beyond)+60]); grid on; set(gca, 'FontSize', 12);
legend('max\_bounces=1', 'max\_bounces=3', 'Location', 'northeast');

subplot(3,1,2);
diff_sig = dp3(ch,:) - dp1(ch,:);
plot(length_domain, diff_sig, 'k-', 'LineWidth', 1.3); hold on;
for j = 1:params.fbg.N_s
    xline(params.fbg.D(j), '-', sprintf('FBG%d (%dG)', j, ghost_counts(j)), ...
        'LineWidth', 1.5, 'Color', colors_g(j,:), 'FontSize', 10, 'LabelOrientation', 'horizontal');
end
for g = 1:length(ghost_beyond)
    xline(ghost_beyond(g), ':', sprintf('G%.0fm', ghost_beyond(g)), ...
        'LineWidth', 1.2, 'Color', [0.8 0 0.6], 'FontSize', 9);
end
xlabel('Pozycja [m]'); ylabel('\Delta Amplituda');
title('Różnica (3b - 1b): ghosty na primary peaks + poza siatkami');
xlim([100, max(ghost_beyond)+60]); grid on; set(gca, 'FontSize', 12);

subplot(3,1,3);
yyaxis left;
b = bar(change_pct, 'FaceColor', 'flat');
b.CData = colors_g;
ylabel('Zmiana amplitudy [%]');
yyaxis right;
plot(1:5, ghost_counts, 'ko-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'k');
ylabel('Liczba ghostów');
xlabel('Numer siatki');
title('Wpływ multi-bounce per siatka');
xticklabels(arrayfun(@(x) sprintf('FBG%d', x), 1:5, 'UniformOutput', false));
grid on; set(gca, 'FontSize', 12);
legend({'Zmiana amp. [%]', 'Liczba ghostów'}, 'Location', 'northwest');

% Zoom na ghosty
figure('Position', [50 50 1200 400], 'theme', 'light', 'Name', 'Zad.1: Zoom ghosty');
plot(length_domain, dp3(ch,:), 'r-', 'LineWidth', 1.5); hold on;
plot(length_domain, dp1(ch,:), 'b--', 'LineWidth', 1);
for g = 1:length(ghost_beyond)
    xline(ghost_beyond(g), ':', sprintf('G%.0fm', ghost_beyond(g)), ...
        'LineWidth', 1.5, 'Color', [0.8 0 0.6], 'FontSize', 10);
end
xlabel('Pozycja [m]'); ylabel('Amplituda');
title('Zoom: ghosty poza ostatnią siatką');
xlim([min(ghost_beyond)-50, max(ghost_beyond)+50]);
grid on; set(gca, 'FontSize', 12);
legend('max\_bounces=3', 'max\_bounces=1', 'Location', 'northeast');

fprintf('\nWykresy gotowe.\n');
