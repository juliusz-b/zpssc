%% Zadanie4_GaussOptimization.m
% Optymalizacja z Gaussian fit - porównanie z centroid
% GA + PSO, 2 warianty (N_s=5, N_s=10)
%
% Projekt: PN/01/0321/2022
% Data: 2025-12-15

clear all; close all; clc;
cd('C:\Users\jbojarczuk\OneDrive - Politechnika Warszawska\PRACA\Insytut Telekomunikacji\Perly Nauki 2022\zpssc');
AddAllSubfolders();
SILENT_MODE = true;
neff = 1.447;

results_dir = fullfile('results', 'WP3', 'optimization');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end

%% ===== KONFIGURACJA =====
p_map = [4, 6, 8, 10];
Nb_map = [2, 4, 8, 16];
nlam_map = [8, 16, 32];
lb = [0.1e-4, 1, 1, 1, 0];
ub = [1e-4, 4, 4, 3, 0.8];
intcon_ga = [2, 3, 4];
nvars = 5;
ga_pop = 10; ga_gen = 8;
pso_swarm = 10; pso_iter = 8;

variants = {
    5,  0,   'N_s=5, grad=0nm';
    10, 0,   'N_s=10, grad=0nm';
};

fprintf('============================================\n');
fprintf('  OPTYMALIZACJA: GAUSSIAN FIT vs CENTROID\n');
fprintf('  %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf('============================================\n\n');

%% ===== BASELINE: Centroid vs Gaussian (default params) =====
fprintf('=== BASELINE ===\n');
for v = 1:size(variants, 1)
    N_s = variants{v,1}; grad_nm = variants{v,2};
    params = defaultParams();
    params.show_plots = false;
    params.fbg.N_s = N_s;
    params.fbg.D = 400:20:(400+20*(N_s-1));
    params.fbg.periods = 530.73e-9*ones(1,N_s);
    params.fbg.deltaneffs = 0.3e-4*ones(1,N_s);
    params.fbg.grating_lengths = 3e-3*ones(1,N_s);
    params.fbg.lambda_shifts = linspace(0, grad_nm*1e-9, N_s);

    out_c = runSimulation(params);
    out_g = runSimulationGauss(params);
    fprintf('  %s: Centroid=%.1f pm, Gauss=%.1f pm\n', variants{v,3}, out_c.MAE, out_g.MAE);
end

%% ===== OPTYMALIZACJA =====
all_gauss = struct();

for v = 1:size(variants, 1)
    N_s = variants{v,1}; grad_nm = variants{v,2}; label = variants{v,3};
    fprintf('\n====== %s (Gaussian fit) ======\n', label);

    objFun = @(x) optimObjectiveGauss(x, p_map, Nb_map, nlam_map, N_s, grad_nm);

    % --- PSO (zazwyczaj lepszy) ---
    fprintf('  PSO (swarm=%d, iter=%d)...\n', pso_swarm, pso_iter);
    opts_pso = optimoptions('particleswarm', ...
        'SwarmSize', pso_swarm, 'MaxIterations', pso_iter, 'Display', 'final');
    tic;
    [x_pso, mae_pso, ~, out_pso] = particleswarm(objFun, nvars, lb, ub, opts_pso);
    t_pso = toc;

    p_pso = p_map(round(x_pso(2))); Nb_pso = Nb_map(round(x_pso(3)));
    nlam_pso = nlam_map(round(x_pso(4)));
    R_pso = tanh(pi*x_pso(1)/(2*neff*530.73e-9)*3e-3)^2*100;
    f_pso = 4e9/((2^p_pso-1)*Nb_pso*nlam_pso*2);

    fprintf('  PSO: MAE=%.1f pm, p=%d, Nb=%d, nlam=%d, dn=%.2e (R=%.1f%%), spread=%.3f\n', ...
        mae_pso, p_pso, Nb_pso, nlam_pso, x_pso(1), R_pso, x_pso(5));
    fprintf('  f_interr=%.1f kHz, %d evals, %.1f s\n', f_pso/1000, out_pso.funccount, t_pso);

    all_gauss(v).label = label;
    all_gauss(v).N_s = N_s;
    all_gauss(v).pso.x = x_pso;
    all_gauss(v).pso.mae = mae_pso;
    all_gauss(v).pso.p = p_pso;
    all_gauss(v).pso.Nb = Nb_pso;
    all_gauss(v).pso.nlam = nlam_pso;
    all_gauss(v).pso.R = R_pso;
    all_gauss(v).pso.f_interr = f_pso;
    all_gauss(v).pso.t = t_pso;
end

%% ===== PORÓWNANIE: Centroid vs Gaussian optimum =====
fprintf('\n============================================\n');
fprintf('  PORÓWNANIE: CENTROID vs GAUSSIAN FIT\n');
fprintf('============================================\n\n');

% Wczytaj wyniki centroid
load(fullfile(results_dir, 'full_optimization_results.mat'), 'all_results');

fprintf('%-25s | %-15s | %-15s | %-15s\n', 'Wariant', 'Default(C)', 'Opt Centroid', 'Opt Gauss');
fprintf('%s\n', repmat('-', 1, 75));
for v = 1:length(all_gauss)
    % Znajdź matching centroid variant
    for vc = 1:length(all_results)
        if all_results(vc).N_s == all_gauss(v).N_s && all_results(vc).grad_nm == 0
            best_c = min(all_results(vc).ga.mae, all_results(vc).pso.mae);
            break;
        end
    end
    def_c = all_results(vc).default_mae;
    gauss_opt = all_gauss(v).pso.mae;
    fprintf('%-25s | %-15.1f | %-15.1f | %-15.1f\n', all_gauss(v).label, def_c, best_c, gauss_opt);
end

%% ===== ZAPIS =====
save(fullfile(results_dir, 'gauss_optimization_results.mat'), 'all_gauss', 'p_map', 'Nb_map', 'nlam_map');

%% ===== WYKRESY =====
figure('Position', [50 50 1000 500], 'theme', 'light', 'Name', 'Centroid vs Gaussian optimization');

mae_compare = zeros(length(all_gauss), 3);
labels_v = {};
for v = 1:length(all_gauss)
    for vc = 1:length(all_results)
        if all_results(vc).N_s == all_gauss(v).N_s && all_results(vc).grad_nm == 0
            mae_compare(v, 1) = all_results(vc).default_mae;
            mae_compare(v, 2) = min(all_results(vc).ga.mae, all_results(vc).pso.mae);
            break;
        end
    end
    mae_compare(v, 3) = all_gauss(v).pso.mae;
    labels_v{v} = all_gauss(v).label;
end

b = bar(mae_compare);
b(1).FaceColor = [0.7 0.7 0.7];
b(2).FaceColor = [0.2 0.5 0.8];
b(3).FaceColor = [0.8 0.2 0.3];
set(gca, 'XTickLabel', labels_v);
ylabel('MAE [pm]');
title(sprintf('Optymalizacja: Default vs Centroid opt vs Gaussian opt - %s', datestr(now)));
legend('Default (centroid)', 'Optimized (centroid)', 'Optimized (Gaussian fit)', 'Location', 'northeast');
grid on; set(gca, 'FontSize', 12);

savefig(gcf, fullfile(results_dir, 'gauss_vs_centroid_optimization.fig'));
saveas(gcf, fullfile(results_dir, 'gauss_vs_centroid_optimization.png'));

fprintf('\nZapisano w %s\n', results_dir);
fprintf('Skrypt zakończony: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
