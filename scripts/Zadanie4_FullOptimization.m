%% Zadanie4_FullOptimization.m
% Pełna optymalizacja GA vs PSO - 5 zmiennych, 4 warianty
% Uruchom cały skrypt - trwa ~10-15 min
%
% Zmienne decyzyjne:
%   x(1): deltaneff    ∈ [0.1e-4, 1e-4]     ciągła
%   x(2): p_idx        ∈ {1,2,3,4} → p=[4,6,8,10]  integer
%   x(3): Nb_idx       ∈ {1,2,3,4} → Nb=[2,4,8,16] integer
%   x(4): nlam_idx     ∈ {1,2,3} → nlam=[8,16,32]   integer
%   x(5): spec_spread  ∈ [0, 0.8] nm         ciągła
%
% Ograniczenie: f_interrogation >= 1 kHz (penalty w funkcji celu)
%
% Warianty:
%   V1: N_s=5,  grad=0nm
%   V2: N_s=10, grad=0nm
%   V3: N_s=10, grad=0.2nm
%   V4: N_s=20, grad=0nm
%
% Projekt: PN/01/0321/2022
% Data: 2025-12-15

clear all; close all; clc;
cd('C:\Users\jbojarczuk\OneDrive - Politechnika Warszawska\PRACA\Insytut Telekomunikacji\Perly Nauki 2022\zpssc');
AddAllSubfolders();

SILENT_MODE = true; % wyłącz disp w generateValues/addNoise

results_dir = fullfile('results', 'WP3', 'optimization');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end

%% ===== KONFIGURACJA =====
p_map = [4, 6, 8, 10];
Nb_map = [2, 4, 8, 16];
nlam_map = [8, 16, 32];

lb = [0.1e-4, 1, 1, 1, 0];
ub = [1e-4,   4, 4, 3, 0.8];
intcon_ga = [2, 3, 4]; % integer variables for GA
nvars = 5;

% Parametry algorytmów
ga_pop = 15;  ga_gen = 12;
pso_swarm = 15; pso_iter = 12;

% Warianty
variants = {
    5,  0,     'N_s=5, grad=0nm';
    10, 0,     'N_s=10, grad=0nm';
    10, 0.2,   'N_s=10, grad=0.2nm';
    20, 0,     'N_s=20, grad=0nm';
};

neff = 1.447;

%% ===== DEFAULT BASELINE =====
fprintf('============================================\n');
fprintf('  PEŁNA OPTYMALIZACJA GA vs PSO\n');
fprintf('  %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf('============================================\n\n');

fprintf('=== BASELINE (default parametry) ===\n');
defaults = {};
for v = 1:size(variants, 1)
    N_s = variants{v,1}; grad_nm = variants{v,2}; label = variants{v,3};
    params = defaultParams();
    params.show_plots = false;
    params.fbg.N_s = N_s;
    params.fbg.D = 400:20:(400+20*(N_s-1));
    params.fbg.periods = 530.73e-9*ones(1,N_s);
    params.fbg.deltaneffs = 0.3e-4*ones(1,N_s);
    params.fbg.grating_lengths = 3e-3*ones(1,N_s);ie
    params.fbg.lambda_shifts = linspace(0, grad_nm*1e-9, N_s);
    tic; out = runSimulation(params); t = toc;
    defaults{v}.mae = out.MAE;
    defaults{v}.label = label;
    fprintf('  %s: MAE = %.1f pm (%.2fs)\n', label, out.MAE, t);
end

%% ===== OPTYMALIZACJA =====
all_results = struct();

for v = 1:size(variants, 1)
    N_s = variants{v,1};
    grad_nm = variants{v,2};
    label = variants{v,3};

    fprintf('\n====== WARIANT %d: %s ======\n', v, label);

    objFun = @(x) optimObjective(x, p_map, Nb_map, nlam_map, N_s, grad_nm);

    % --- GA ---
    fprintf('\n  GA (pop=%d, gen=%d)...\n', ga_pop, ga_gen);
    opts_ga = optimoptions('ga', ...
        'PopulationSize', ga_pop, ...
        'MaxGenerations', ga_gen, ...
        'Display', 'final', ...
        'UseParallel', false);

    tic;
    [x_ga, mae_ga, ~, out_ga] = ga(objFun, nvars, [], [], [], [], lb, ub, [], intcon_ga, opts_ga);
    t_ga = toc;

    p_ga = p_map(x_ga(2)); Nb_ga = Nb_map(x_ga(3)); nlam_ga = nlam_map(x_ga(4));
    f_ga = 4e9 / ((2^p_ga-1) * Nb_ga * nlam_ga * 2);
    R_ga = tanh(pi*x_ga(1)/(2*neff*530.73e-9)*3e-3)^2*100;

    fprintf('  GA wynik: MAE=%.1f pm, dn=%.2e (R=%.1f%%), p=%d, Nb=%d, nlam=%d, spread=%.3f\n', ...
        mae_ga, x_ga(1), R_ga, p_ga, Nb_ga, nlam_ga, x_ga(5));
    fprintf('  GA: %d ewaluacji, %.1f s, f_interr=%.1f kHz\n', out_ga.funccount, t_ga, f_ga/1000);

    % --- PSO ---
    fprintf('\n  PSO (swarm=%d, iter=%d)...\n', pso_swarm, pso_iter);
    opts_pso = optimoptions('particleswarm', ...
        'SwarmSize', pso_swarm, ...
        'MaxIterations', pso_iter, ...
        'Display', 'final');

    tic;
    [x_pso, mae_pso, ~, out_pso] = particleswarm(objFun, nvars, lb, ub, opts_pso);
    t_pso = toc;

    p_pso = p_map(round(x_pso(2))); Nb_pso = Nb_map(round(x_pso(3))); nlam_pso = nlam_map(round(x_pso(4)));
    f_pso = 4e9 / ((2^p_pso-1) * Nb_pso * nlam_pso * 2);
    R_pso = tanh(pi*x_pso(1)/(2*neff*530.73e-9)*3e-3)^2*100;

    fprintf('  PSO wynik: MAE=%.1f pm, dn=%.2e (R=%.1f%%), p=%d, Nb=%d, nlam=%d, spread=%.3f\n', ...
        mae_pso, x_pso(1), R_pso, p_pso, Nb_pso, nlam_pso, x_pso(5));
    fprintf('  PSO: %d ewaluacji, %.1f s, f_interr=%.1f kHz\n', out_pso.funccount, t_pso, f_pso/1000);

    % Zapisz wyniki wariantu
    all_results(v).label = label;
    all_results(v).N_s = N_s;
    all_results(v).grad_nm = grad_nm;
    all_results(v).default_mae = defaults{v}.mae;
    all_results(v).ga.x = x_ga; all_results(v).ga.mae = mae_ga;
    all_results(v).ga.t = t_ga; all_results(v).ga.funccount = out_ga.funccount;
    all_results(v).ga.p = p_ga; all_results(v).ga.Nb = Nb_ga;
    all_results(v).ga.nlam = nlam_ga; all_results(v).ga.R = R_ga;
    all_results(v).ga.f_interr = f_ga;
    all_results(v).pso.x = x_pso; all_results(v).pso.mae = mae_pso;
    all_results(v).pso.t = t_pso; all_results(v).pso.funccount = out_pso.funccount;
    all_results(v).pso.p = p_pso; all_results(v).pso.Nb = Nb_pso;
    all_results(v).pso.nlam = nlam_pso; all_results(v).pso.R = R_pso;
    all_results(v).pso.f_interr = f_pso;
end

%% ===== ZAPIS =====
save(fullfile(results_dir, 'full_optimization_results.mat'), ...
    'all_results', 'p_map', 'Nb_map', 'nlam_map', 'variants', ...
    'ga_pop', 'ga_gen', 'pso_swarm', 'pso_iter');

%% ===== PODSUMOWANIE =====
fprintf('\n\n============================================\n');
fprintf('  PODSUMOWANIE OPTYMALIZACJI\n');
fprintf('  %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf('============================================\n\n');

fprintf('%-25s | %-12s | %-30s | %-30s\n', 'Wariant', 'Default', 'GA', 'PSO');
fprintf('%s\n', repmat('-', 1, 105));
for v = 1:length(all_results)
    r = all_results(v);
    ga_str = sprintf('MAE=%.0f, p=%d, Nb=%d, sp=%.2f', r.ga.mae, r.ga.p, r.ga.Nb, r.ga.x(5));
    pso_str = sprintf('MAE=%.0f, p=%d, Nb=%d, sp=%.2f', r.pso.mae, r.pso.p, r.pso.Nb, r.pso.x(5));
    fprintf('%-25s | MAE=%-7.0f | %-30s | %-30s\n', r.label, r.default_mae, ga_str, pso_str);
end

% Poprawa procentowa
fprintf('\n%-25s | %-12s | %-12s\n', 'Wariant', 'Poprawa GA', 'Poprawa PSO');
fprintf('%s\n', repmat('-', 1, 55));
for v = 1:length(all_results)
    r = all_results(v);
    imp_ga = (r.default_mae - r.ga.mae) / r.default_mae * 100;
    imp_pso = (r.default_mae - r.pso.mae) / r.default_mae * 100;
    fprintf('%-25s | %+.1f%%        | %+.1f%%\n', r.label, imp_ga, imp_pso);
end

%% ===== WYKRESY =====

% Wykres 1: Bar chart - Default vs GA vs PSO per wariant
figure('Position', [50 50 1200 500], 'theme', 'light', 'Name', 'Zad.4: Optymalizacja GA vs PSO');

mae_data = zeros(length(all_results), 3);
labels_v = {};
for v = 1:length(all_results)
    mae_data(v, :) = [all_results(v).default_mae, all_results(v).ga.mae, all_results(v).pso.mae];
    labels_v{v} = all_results(v).label;
end

subplot(1,2,1);
b = bar(mae_data);
b(1).FaceColor = [0.7 0.7 0.7]; % default: szary
b(2).FaceColor = [0.2 0.5 0.8]; % GA: niebieski
b(3).FaceColor = [0.8 0.3 0.2]; % PSO: czerwony
set(gca, 'XTickLabel', labels_v, 'XTickLabelRotation', 20);
ylabel('MAE [pm]'); title('MAE: Default vs GA vs PSO');
legend('Default', 'GA', 'PSO', 'Location', 'northwest');
grid on; set(gca, 'FontSize', 11);

% Wykres 2: Optymalne parametry - spider/radar lub tabela wizualna
subplot(1,2,2);
% Poprawa % per wariant
imp_ga = zeros(1, length(all_results));
imp_pso = zeros(1, length(all_results));
for v = 1:length(all_results)
    imp_ga(v) = (all_results(v).default_mae - all_results(v).ga.mae) / all_results(v).default_mae * 100;
    imp_pso(v) = (all_results(v).default_mae - all_results(v).pso.mae) / all_results(v).default_mae * 100;
end
bar_imp = [imp_ga; imp_pso]';
b2 = bar(bar_imp);
b2(1).FaceColor = [0.2 0.5 0.8];
b2(2).FaceColor = [0.8 0.3 0.2];
set(gca, 'XTickLabel', labels_v, 'XTickLabelRotation', 20);
ylabel('Poprawa vs default [%]'); title('Poprawa optymalizacji');
legend('GA', 'PSO', 'Location', 'northwest');
grid on; set(gca, 'FontSize', 11);

sgtitle(sprintf('Zadanie 4: Optymalizacja GA vs PSO - %s', datestr(now, 'yyyy-mm-dd HH:MM')), ...
    'FontSize', 14, 'FontWeight', 'bold');

savefig(gcf, fullfile(results_dir, 'full_optimization_comparison.fig'));
saveas(gcf, fullfile(results_dir, 'full_optimization_comparison.png'));

% Wykres 3: Optymalne parametry per wariant - szczegóły
figure('Position', [50 600 1200 400], 'theme', 'light', 'Name', 'Zad.4: Optymalne parametry');
param_names = {'deltaneff [×1e-4]', 'p', 'Nb', 'n\_lambdas', 'spread [nm]'};
for v = 1:length(all_results)
    r = all_results(v);
    subplot(1, length(all_results), v);
    vals_def = [0.3, 8, 8, 16, 0];
    vals_ga = [r.ga.x(1)*1e4, r.ga.p, r.ga.Nb, r.ga.nlam, r.ga.x(5)];
    vals_pso = [r.pso.x(1)*1e4, r.pso.p, r.pso.Nb, r.pso.nlam, r.pso.x(5)];

    data_plot = [vals_def; vals_ga; vals_pso]';
    barh(data_plot);
    set(gca, 'YTickLabel', param_names);
    xlabel('Wartość'); title(sprintf('%s\nGA=%.0f, PSO=%.0f pm', r.label, r.ga.mae, r.pso.mae));
    legend('Default', 'GA', 'PSO', 'Location', 'best', 'FontSize', 8);
    set(gca, 'FontSize', 10);
end

savefig(gcf, fullfile(results_dir, 'full_optimization_params.fig'));
saveas(gcf, fullfile(results_dir, 'full_optimization_params.png'));

fprintf('\nWszystko zapisane w %s\n', results_dir);
fprintf('Skrypt zakończony: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
