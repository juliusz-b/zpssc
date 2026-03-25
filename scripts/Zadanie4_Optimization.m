%% Zadanie4_Optimization.m
% Optymalizacja parametrów systemu CDM-FBG (WP3.2 + WP3.3)
% Grid search: p × deltaneff dla różnych wariantów (N_s, gradient temp.)
%
% Zmienne decyzyjne:
%   p ∈ {4, 6, 8, 10} (Kasami, parzyste)
%   deltaneff ∈ [0.1e-4, 1e-4] (8 punktów)
%
% Funkcja celu: min MAE (mean |Δλ_B|) [pm]
%
% Projekt: PN/01/0321/2022
% Data: 2025-12-15

clear all; close all; clc;
cd('C:\Users\jbojarczuk\OneDrive - Politechnika Warszawska\PRACA\Insytut Telekomunikacji\Perly Nauki 2022\zpssc');
AddAllSubfolders();

results_dir = fullfile('results', 'WP3', 'optimization');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end

%% ===== PARAMETRY GRID =====
p_vals = [4, 6, 8, 10];
dn_vals = linspace(0.1e-4, 1e-4, 8);
neff = 1.447;

kappa_fn = @(dn) pi*dn./(2*neff*530.73e-9);
Rpk_fn = @(dn) tanh(kappa_fn(dn)*3e-3).^2 * 100;

%% ===== WARIANTY =====
variants = {
    10, 0,    'N_s=10, grad=0nm';
    10, 0.2,  'N_s=10, grad=0.2nm';
};

all_results = {};
fprintf('=== OPTYMALIZACJA: %d wariantów × %d grid points ===\n\n', ...
    size(variants,1), length(p_vals)*length(dn_vals));

for v = 1:size(variants, 1)
    N_s = variants{v, 1};
    grad_nm = variants{v, 2};
    label = variants{v, 3};

    fprintf('--- %s ---\n', label);
    mae_grid = NaN(length(p_vals), length(dn_vals));
    tic;

    for ip = 1:length(p_vals)
        for idn = 1:length(dn_vals)
            params = defaultParams();
            params.show_plots = false;
            params.p = p_vals(ip);
            params.fbg.N_s = N_s;
            params.fbg.D = 400:20:(400+20*(N_s-1));
            params.fbg.periods = 530.73e-9*ones(1,N_s);
            params.fbg.deltaneffs = dn_vals(idn)*ones(1,N_s);
            params.fbg.grating_lengths = 3e-3*ones(1,N_s);
            params.fbg.lambda_shifts = linspace(0, grad_nm*1e-9, N_s);
            out = runSimulation(params);
            mae_grid(ip, idn) = out.MAE;
        end
    end

    [min_mae, idx] = min(mae_grid(:));
    [ip_opt, idn_opt] = ind2sub(size(mae_grid), idx);

    all_results{v}.mae_grid = mae_grid;
    all_results{v}.p_opt = p_vals(ip_opt);
    all_results{v}.dn_opt = dn_vals(idn_opt);
    all_results{v}.mae_opt = min_mae;
    all_results{v}.N_s = N_s;
    all_results{v}.grad_nm = grad_nm;
    all_results{v}.label = label;

    fprintf('  Optimum: p=%d, deltaneff=%.2e (R=%.1f%%), MAE=%.1f pm\n', ...
        p_vals(ip_opt), dn_vals(idn_opt), Rpk_fn(dn_vals(idn_opt)), min_mae);
    fprintf('  Czas: %.1f s\n\n', toc);
end

%% ===== ZAPIS =====
save(fullfile(results_dir, 'optimization_results.mat'), ...
    'all_results', 'p_vals', 'dn_vals', 'variants');

%% ===== WYKRESY =====
rpk_vals = Rpk_fn(dn_vals);

figure('Position', [50 50 1300 600], 'theme', 'light', 'Name', 'Zad.4: Optymalizacja');
for v = 1:size(variants, 1)
    subplot(1, size(variants,1), v);
    imagesc(rpk_vals, p_vals, all_results{v}.mae_grid);
    set(gca, 'YDir', 'normal');
    cb = colorbar; cb.Label.String = 'MAE [pm]';
    colormap(flipud(hot));
    xlabel('R_{peak} [%]'); ylabel('p');
    title(sprintf('(%c) %s\nOpt: p=%d, R=%.1f%%, MAE=%.0f pm', ...
        char('a'+v-1), all_results{v}.label, all_results{v}.p_opt, ...
        Rpk_fn(all_results{v}.dn_opt), all_results{v}.mae_opt));
    yticks(p_vals); set(gca, 'FontSize', 12);
    hold on;
    [~, ri] = min(abs(dn_vals - all_results{v}.dn_opt));
    plot(rpk_vals(ri), all_results{v}.p_opt, 'gp', 'MarkerSize', 20, ...
        'MarkerFaceColor', 'g', 'LineWidth', 2);
end
sgtitle('Optymalizacja p × \Deltan_{eff}: heatmap MAE [pm]', 'FontSize', 14, 'FontWeight', 'bold');

savefig(gcf, fullfile(results_dir, 'optimization_heatmap.fig'));
saveas(gcf, fullfile(results_dir, 'optimization_heatmap.png'));

%% ===== PODSUMOWANIE =====
fprintf('\n========================================\n');
fprintf('  WYNIKI OPTYMALIZACJI (2025-12-15)\n');
fprintf('========================================\n');
fprintf('%-25s %-6s %-12s %-10s %-10s\n', 'Wariant', 'p_opt', 'dn_opt', 'R[%%]', 'MAE[pm]');
for v = 1:length(all_results)
    fprintf('%-25s %-6d %-12.2e %-10.1f %-10.1f\n', ...
        all_results{v}.label, all_results{v}.p_opt, all_results{v}.dn_opt, ...
        Rpk_fn(all_results{v}.dn_opt), all_results{v}.mae_opt);
end
fprintf('\nWszystko zapisane w %s\n', results_dir);
