%% Zadanie3_SensitivityAnalysis.m
% Analiza wrażliwości (WP3.1) - skanowanie parametrów po jednym na raz
%
% Projekt: PN/01/0321/2022

clear all; close all; clc;
cd('C:\Users\jbojarczuk\OneDrive - Politechnika Warszawska\PRACA\Insytut Telekomunikacji\Perly Nauki 2022\zpssc');
AddAllSubfolders();

results_dir = fullfile('results', 'WP3', 'sensitivity');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end

neff = 1.447;

%% ===== SWEEP 1: p (długość kodu) =====
fprintf('=== SWEEP: p ===\n');
p_vals = [4, 6, 8, 10]; mae_p = zeros(size(p_vals));
for idx = 1:length(p_vals)
    params = defaultParams();
    params.show_plots = false;
    params.p = p_vals(idx);
    fprintf('p=%d... ', p_vals(idx)); tic;
    out = runSimulation(params);
    mae_p(idx) = out.MAE;
    fprintf('MAE=%.1f pm, t=%.1fs\n', out.MAE, toc);
end

%% ===== SWEEP 2: deltaneff =====
fprintf('\n=== SWEEP: deltaneff ===\n');
dn_vals = [0.1e-4, 0.3e-4, 0.5e-4, 0.7e-4, 1e-4];
mae_dn = zeros(size(dn_vals));
rpk = zeros(size(dn_vals));
for idx = 1:length(dn_vals)
    dn = dn_vals(idx);
    params = defaultParams();
    params.show_plots = false;
    params.fbg.deltaneffs = dn * ones(1, params.fbg.N_s);
    kappa = pi*dn/(2*neff*params.fbg.periods(1));
    rpk(idx) = tanh(kappa*3e-3)^2 * 100;
    fprintf('dneff=%.1e (R=%.1f%%)... ', dn, rpk(idx)); tic;
    out = runSimulation(params);
    mae_dn(idx) = out.MAE;
    fprintf('MAE=%.1f pm, t=%.1fs\n', out.MAE, toc);
end

%% ===== SWEEP 3: N_s =====
fprintf('\n=== SWEEP: N_s ===\n');
ns_vals = [3, 5, 10, 15]; mae_ns = zeros(size(ns_vals));
for idx = 1:length(ns_vals)
    ns = ns_vals(idx);
    params = defaultParams();
    params.show_plots = false;
    params.fbg.N_s = ns;
    params.fbg.D = 400:20:(400+20*(ns-1));
    params.fbg.periods = 530.73e-9 * ones(1, ns);
    params.fbg.deltaneffs = 0.3e-4 * ones(1, ns);
    params.fbg.grating_lengths = 3e-3 * ones(1, ns);
    params.fbg.lambda_shifts = zeros(1, ns);
    fprintf('N_s=%d... ', ns); tic;
    out = runSimulation(params);
    mae_ns(idx) = out.MAE;
    fprintf('MAE=%.1f pm, t=%.1fs\n', out.MAE, toc);
end

%% ===== SWEEP 4: NEP (słaby sygnał: N_s=10, deltaneff=0.1e-4) =====
fprintf('\n=== SWEEP: NEP (N_s=10, R=0.38%%) ===\n');
nep_vals = [1, 50, 200]; mae_nep = zeros(size(nep_vals));
for idx = 1:length(nep_vals)
    params = defaultParams();
    params.show_plots = false;
    params.fbg.N_s = 10;
    params.fbg.D = 400:20:(400+20*9);
    params.fbg.periods = 530.73e-9 * ones(1, 10);
    params.fbg.deltaneffs = 0.1e-4 * ones(1, 10);
    params.fbg.grating_lengths = 3e-3 * ones(1, 10);
    params.fbg.lambda_shifts = zeros(1, 10);
    params.noise.NEP = nep_vals(idx) * 1e-12;
    fprintf('NEP=%dpW... ', nep_vals(idx)); tic;
    out = runSimulation(params);
    mae_nep(idx) = out.MAE;
    fprintf('MAE=%.1f pm, t=%.1fs\n', out.MAE, toc);
end

%% ===== SWEEP 5: gradient temperaturowy =====
fprintf('\n=== SWEEP: gradient temperaturowy ===\n');
gr_vals = [0, 0.1, 0.2, 0.5]; mae_gr = zeros(size(gr_vals));
for idx = 1:length(gr_vals)
    params = defaultParams();
    params.show_plots = false;
    params.fbg.lambda_shifts = linspace(0, gr_vals(idx)*1e-9, params.fbg.N_s);
    fprintf('grad=%.1fnm... ', gr_vals(idx)); tic;
    out = runSimulation(params);
    mae_gr(idx) = out.MAE;
    fprintf('MAE=%.1f pm, t=%.1fs\n', out.MAE, toc);
end

%% ===== ZAPIS =====
save(fullfile(results_dir, 'sensitivity_results.mat'), ...
    'p_vals', 'mae_p', 'dn_vals', 'mae_dn', 'rpk', 'ns_vals', 'mae_ns', ...
    'gr_vals', 'mae_gr', 'nep_vals', 'mae_nep');
fprintf('\nDane zapisane w %s\n', results_dir);

%% ===== WYKRESY =====
figure('Position', [50 50 1400 900], 'theme', 'light', 'Name', 'Zad.3: Analiza wrażliwości (WP3.1)');

subplot(2,3,1);
plot(p_vals, mae_p, 'bo-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
xlabel('p'); ylabel('MAE [pm]');
title('(a) Długość kodu p'); grid on; set(gca, 'FontSize', 11);

subplot(2,3,2);
plot(dn_vals*1e4, mae_dn, 'ro-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('\Deltan_{eff} [\times10^{-4}]'); ylabel('MAE [pm]');
title('(b) Głębokość modulacji'); grid on; set(gca, 'FontSize', 11);

subplot(2,3,3);
plot(ns_vals, mae_ns, 'gs-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0 0.7 0]);
xlabel('N_s'); ylabel('MAE [pm]');
title('(c) Liczba siatek'); grid on; set(gca, 'FontSize', 11);

subplot(2,3,4);
plot(gr_vals, mae_gr, 'md-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'm');
xlabel('Gradient [nm]'); ylabel('MAE [pm]');
title('(d) Gradient temperaturowy'); grid on; set(gca, 'FontSize', 11);

subplot(2,3,5);
plot(nep_vals, mae_nep, 'ks-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0.5 0.5 0.5]);
xlabel('NEP [pW/\surdHz]'); ylabel('MAE [pm]');
title('(e) NEP (N_s=10, R=0.38%)'); grid on; set(gca, 'FontSize', 11);

subplot(2,3,6);
ranges = [max(mae_p)-min(mae_p); max(mae_dn)-min(mae_dn); max(mae_ns)-min(mae_ns);
          max(mae_gr)-min(mae_gr); max(mae_nep)-min(mae_nep)];
param_names = {'p (kod)', '\Deltan_{eff}', 'N_s', 'Gradient T', 'NEP'};
[ranges_s, si] = sort(ranges, 'descend');
bar(ranges_s, 'FaceColor', [0.3 0.5 0.8]);
set(gca, 'XTickLabel', param_names(si), 'XTickLabelRotation', 20);
ylabel('\DeltaMAE [pm]'); title('(f) Ranking wrażliwości');
grid on; set(gca, 'FontSize', 11);

sgtitle('Zadanie 3: Analiza wrażliwości (WP3.1)', 'FontSize', 14, 'FontWeight', 'bold');

savefig(gcf, fullfile(results_dir, 'sensitivity_all.fig'));
saveas(gcf, fullfile(results_dir, 'sensitivity_all.png'));
fprintf('Wykresy gotowe.\n');
