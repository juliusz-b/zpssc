%% Zadanie_PINvsAPD.m
% Porównanie fotodetektora PIN vs APD w systemie CDM-FBG
% Analiza wpływu gainu APD na MAE detekcji λ_B
%
% Model APD: excess noise factor F(M) = k_eff*M + (1-k_eff)*(2-1/M)
% InGaAs: k_eff ≈ 0.7
% Efektywne NEP: NEP_eff = NEP_base * sqrt(F(M)) / M
%
% Projekt: PN/01/0321/2022

clear all; close all; clc;
cd('C:\Users\jbojarczuk\OneDrive - Politechnika Warszawska\PRACA\Insytut Telekomunikacji\Perly Nauki 2022\zpssc');
AddAllSubfolders();

results_dir = fullfile('results', 'WP3', 'sensitivity');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end

%% ===== PARAMETRY APD =====
k_eff = 0.7;  % współczynnik jonizacji InGaAs
F = @(M) k_eff * M + (1 - k_eff) * (2 - 1./M);
gains = [1, 2, 5, 10, 20, 50];  % M=1 = PIN

%% ===== SCENARIUSZ 1: Typowe warunki =====
NEP_base1 = 15e-12;
fprintf('=== Scenariusz 1: Typowe warunki (N_s=5, R=3.3%%, NEP=15pW) ===\n');
mae1 = zeros(size(gains));
nep1 = zeros(size(gains));
for idx = 1:length(gains)
    M = gains(idx);
    nep1(idx) = NEP_base1 * sqrt(F(M)) / M;
    params = defaultParams();
    params.show_plots = false;
    params.pd.gain = M;
    params.noise.NEP = nep1(idx);
    out = runSimulation(params);
    mae1(idx) = out.MAE;
    fprintf('  M=%-3d: NEP_eff=%.1f pW, MAE=%.1f pm\n', M, nep1(idx)*1e12, mae1(idx));
end

%% ===== SCENARIUSZ 2: Warunki niekorzystne =====
NEP_base2 = 50e-12;
fprintf('\n=== Scenariusz 2: Niekorzystne (N_s=10, R=0.38%%, NEP=50pW) ===\n');
mae2 = zeros(size(gains));
nep2 = zeros(size(gains));
for idx = 1:length(gains)
    M = gains(idx);
    nep2(idx) = NEP_base2 * sqrt(F(M)) / M;
    params = defaultParams();
    params.show_plots = false;
    params.pd.gain = M;
    params.noise.NEP = nep2(idx);
    params.fbg.N_s = 10;
    params.fbg.D = 400:20:(400+20*9);
    params.fbg.periods = 530.73e-9 * ones(1, 10);
    params.fbg.deltaneffs = 0.1e-4 * ones(1, 10);
    params.fbg.grating_lengths = 3e-3 * ones(1, 10);
    params.fbg.lambda_shifts = zeros(1, 10);
    out = runSimulation(params);
    mae2(idx) = out.MAE;
    fprintf('  M=%-3d: NEP_eff=%.1f pW, MAE=%.1f pm\n', M, nep2(idx)*1e12, mae2(idx));
end

%% ===== WYNIKI =====
fprintf('\n=== PODSUMOWANIE ===\n');
fprintf('Scenariusz 1 (typowy): PIN MAE=%.1f pm, APD(M=10) MAE=%.1f pm → zmiana %+.1f pm\n', ...
    mae1(1), mae1(gains==10), mae1(gains==10)-mae1(1));
fprintf('Scenariusz 2 (niekorzystny): PIN MAE=%.1f pm, APD(M=10) MAE=%.1f pm → zmiana %+.1f pm (%.1f%%)\n', ...
    mae2(1), mae2(gains==10), mae2(gains==10)-mae2(1), (mae2(1)-mae2(gains==10))/mae2(1)*100);

%% ===== WYKRESY =====
figure('Position', [50 50 1200 500], 'theme', 'light', 'Name', 'PIN vs APD');

subplot(1,2,1);
plot(gains, mae1, 'bo-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'b'); hold on;
plot(gains, mae2, 'rs-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
set(gca, 'XScale', 'log');
xlabel('APD Gain M'); ylabel('MAE [pm]');
title('PIN vs APD: MAE vs gain');
legend('Typowe (N_s=5, R=3.3%, NEP=15pW)', 'Niekorzystne (N_s=10, R=0.4%, NEP=50pW)', ...
    'Location', 'best');
grid on; set(gca, 'FontSize', 11);
xticks(gains); xticklabels(arrayfun(@(x) sprintf('M=%d', x), gains, 'UniformOutput', false));

subplot(1,2,2);
plot(gains, nep1*1e12, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b'); hold on;
plot(gains, nep2*1e12, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
set(gca, 'XScale', 'log');
xlabel('APD Gain M'); ylabel('NEP_{eff} [pW/\surdHz]');
title(sprintf('Efektywne NEP vs gain (k_{eff}=%.1f)', k_eff));
legend('NEP_{base}=15 pW', 'NEP_{base}=50 pW', 'Location', 'best');
grid on; set(gca, 'FontSize', 11);
xticks(gains); xticklabels(arrayfun(@(x) sprintf('M=%d', x), gains, 'UniformOutput', false));

savefig(gcf, fullfile(results_dir, 'pin_vs_apd.fig'));
saveas(gcf, fullfile(results_dir, 'pin_vs_apd.png'));

% Zapis danych
save(fullfile(results_dir, 'pin_vs_apd_results.mat'), ...
    'gains', 'mae1', 'mae2', 'nep1', 'nep2', 'k_eff', 'NEP_base1', 'NEP_base2');

%% ===== SCENARIUSZ 3: MAE vs moc optyczna (parametry Thorlabs) =====
% PIN: R=0.9 A/W, NEP=15 pW/√Hz, M=1 (np. Thorlabs PDA10CS2)
% APD: R=0.9 A/W base, M=10, NEP=2 pW/√Hz input-referred (np. Thorlabs APD430C)

pin_th.R = 0.9; pin_th.NEP = 15e-12; pin_th.gain = 1;
apd_th.R = 0.9; apd_th.NEP = 2e-12;  apd_th.gain = 10;

P_dBm = [-30, -25, -20, -18, -15, -12, -10, -5, 0];
P_W = 10.^(P_dBm/10) * 1e-3;
mae_pin_p = zeros(size(P_dBm));
mae_apd_p = zeros(size(P_dBm));

fprintf('\n=== MAE vs moc optyczna (Thorlabs) ===\n');
for idx = 1:length(P_dBm)
    % PIN
    params = defaultParams();
    params.show_plots = false;
    params.laser.power = P_W(idx);
    params.pd.A = pin_th.R; params.pd.gain = pin_th.gain; params.noise.NEP = pin_th.NEP;
    out = runSimulation(params);
    mae_pin_p(idx) = out.MAE;

    % APD
    params.pd.A = apd_th.R; params.pd.gain = apd_th.gain; params.noise.NEP = apd_th.NEP;
    out = runSimulation(params);
    mae_apd_p(idx) = out.MAE;

    fprintf('  P=%+3d dBm (%6.1f µW): PIN=%.1f pm, APD=%.1f pm\n', ...
        P_dBm(idx), P_W(idx)*1e6, mae_pin_p(idx), mae_apd_p(idx));
end

% Wykres: MAE vs moc optyczna
cross_idx = find(mae_pin_p < mae_apd_p, 1, 'first');
cross_dBm = P_dBm(cross_idx);

figure('Position', [50 50 900 600], 'theme', 'light', 'Name', 'PIN vs APD: MAE vs moc optyczna');
plot(P_dBm, mae_pin_p, 'b^-', 'LineWidth', 2.5, 'MarkerSize', 12, 'MarkerFaceColor', 'b'); hold on;
plot(P_dBm, mae_apd_p, 'rv-', 'LineWidth', 2.5, 'MarkerSize', 12, 'MarkerFaceColor', 'r');
xline(cross_dBm, '--k', sprintf('Crossing ≈ %d dBm', cross_dBm), ...
    'LineWidth', 1.5, 'FontSize', 12, 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
fill([-32 cross_dBm cross_dBm -32], [330 330 470 470], [1 0.9 0.9], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
fill([cross_dBm 2 2 cross_dBm], [330 330 470 470], [0.9 0.9 1], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
text(-27, 460, 'APD lepszy', 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.8 0 0], 'HorizontalAlignment', 'center');
text(-5, 460, 'PIN lepszy (lub równy)', 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0 0 0.7], 'HorizontalAlignment', 'center');
yline(393, ':k', 'Floor sidelobes', 'LineWidth', 1, 'FontSize', 10, 'LabelHorizontalAlignment', 'left');
xlabel('Moc optyczna lasera [dBm]', 'FontSize', 14);
ylabel('MAE [pm]', 'FontSize', 14);
title(sprintf('PIN (NEP=15pW) vs APD Thorlabs (M=10, NEP=2pW) - %d FBG, p=%d', 5, 8), 'FontSize', 14);
legend('PIN (R=0.9 A/W, NEP=15 pW/\surdHz)', 'APD (M=10, NEP=2 pW/\surdHz)', 'Location', 'northeast', 'FontSize', 11);
grid on; set(gca, 'FontSize', 12);
xlim([-32, 2]); ylim([330, 470]);

savefig(gcf, fullfile(results_dir, 'pin_vs_apd_power.fig'));
saveas(gcf, fullfile(results_dir, 'pin_vs_apd_power.png'));

% Zapis
save(fullfile(results_dir, 'pin_vs_apd_power_results.mat'), ...
    'P_dBm', 'P_W', 'mae_pin_p', 'mae_apd_p', 'pin_th', 'apd_th');

fprintf('\nWszystko zapisane w %s\n', results_dir);
