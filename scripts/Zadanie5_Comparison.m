%% Zadanie5_Comparison.m
% Porównanie przed/po optymalizacji (WP3.4)
% Wczytuje wyniki z Zadanie4_FullOptimization.m i generuje wykresy porównawcze
%
% Wymaga: results/WP3/optimization/full_optimization_results.mat
%         (wygenerowane przez Zadanie4_FullOptimization.m)
%
% Projekt: PN/01/0321/2022
% Data: 2025-12-15

clear all; close all; clc;
cd('C:\Users\jbojarczuk\OneDrive - Politechnika Warszawska\PRACA\Insytut Telekomunikacji\Perly Nauki 2022\zpssc');
AddAllSubfolders();

SILENT_MODE = true;
neff = 1.447;

results_dir = fullfile('results', 'WP3', 'comparison');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end

%% ===== WCZYTAJ WYNIKI OPTYMALIZACJI =====
opt_file = fullfile('results', 'WP3', 'optimization', 'full_optimization_results.mat');
if ~exist(opt_file, 'file')
    error('Brak pliku %s - najpierw uruchom Zadanie4_FullOptimization.m', opt_file);
end
load(opt_file);

fprintf('============================================\n');
fprintf('  PORÓWNANIE PRZED/PO OPTYMALIZACJI (WP3.4)\n');
fprintf('  %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf('============================================\n\n');

%% ===== PORÓWNANIE DLA KAŻDEGO WARIANTU =====
% Dla najlepszego algorytmu (PSO lub GA - ten z niższym MAE) uruchom:
% 1. Symulację default
% 2. Symulację z optymalnymi parametrami
% 3. Porównaj: korelacje, widma, błędy per siatka

for v = 1:length(all_results)
    r = all_results(v);
    N_s = r.N_s;
    grad_nm = r.grad_nm;

    fprintf('\n====== WARIANT %d: %s ======\n', v, r.label);

    % Wybierz lepszy algorytm
    if r.pso.mae <= r.ga.mae
        best = r.pso;
        best_name = 'PSO';
    else
        best = r.ga;
        best_name = 'GA';
    end

    fprintf('Lepszy algorytm: %s (MAE=%.1f pm vs %.1f pm)\n', best_name, best.mae, ...
        ternary(strcmp(best_name,'PSO'), r.ga.mae, r.pso.mae));

    % --- Symulacja DEFAULT ---
    params_def = defaultParams();
    params_def.show_plots = false;
    params_def.fbg.N_s = N_s;
    params_def.fbg.D = 400:20:(400+20*(N_s-1));
    params_def.fbg.periods = 530.73e-9*ones(1,N_s);
    params_def.fbg.deltaneffs = 0.3e-4*ones(1,N_s);
    params_def.fbg.grating_lengths = 3e-3*ones(1,N_s);
    params_def.fbg.lambda_shifts = linspace(0, grad_nm*1e-9, N_s);

    out_def = runSimulation(params_def);

    % --- Symulacja OPTIMUM ---
    params_opt = defaultParams();
    params_opt.show_plots = false;
    params_opt.p = best.p;
    params_opt.Nb = best.Nb;
    params_opt.laser.step = best.nlam;
    params_opt.modulation.Fn_s = (20:90:(20+90*(best.nlam-1)))*1e6;
    params_opt.fbg.N_s = N_s;
    params_opt.fbg.D = 400:20:(400+20*(N_s-1));
    params_opt.fbg.periods = 530.73e-9*ones(1,N_s);
    params_opt.fbg.deltaneffs = best.x(1)*ones(1,N_s);
    params_opt.fbg.grating_lengths = 3e-3*ones(1,N_s);
    % Spectral spread + gradient
    params_opt.fbg.lambda_shifts = linspace(-best.x(5)/2, best.x(5)/2, N_s)*1e-9 ...
        + linspace(0, grad_nm*1e-9, N_s);

    out_opt = runSimulation(params_opt);

    % --- Wyniki tekstowe ---
    fprintf('\n  %-20s %-15s %-15s\n', '', 'DEFAULT', best_name);
    fprintf('  %-20s %-15s %-15s\n', '---', '-------', '-------');
    fprintf('  %-20s %-15d %-15d\n', 'p', 8, best.p);
    fprintf('  %-20s %-15d %-15d\n', 'Nb', 8, best.Nb);
    fprintf('  %-20s %-15d %-15d\n', 'n_lambdas', 16, best.nlam);
    fprintf('  %-20s %-15.2e %-15.2e\n', 'deltaneff', 0.3e-4, best.x(1));
    fprintf('  %-20s %-15.1f %-15.1f\n', 'R_peak [%%]', ...
        tanh(pi*0.3e-4/(2*neff*530.73e-9)*3e-3)^2*100, best.R);
    fprintf('  %-20s %-15.3f %-15.3f\n', 'spec_spread [nm]', 0, best.x(5));
    fprintf('  %-20s %-15.1f %-15.1f\n', 'f_interr [kHz]', ...
        4e9/((2^8-1)*8*16*2)/1000, best.f_interr/1000);
    fprintf('  %-20s %-15.1f %-15.1f\n', 'MAE [pm]', out_def.MAE, out_opt.MAE);
    fprintf('  %-20s %-15s %-15s\n', 'Poprawa', '', ...
        sprintf('-%.0f%%', (out_def.MAE-out_opt.MAE)/out_def.MAE*100));

    % --- Wykresy ---
    figure('Position', [50 50+200*(v-1) 1400 600], 'theme', 'light', ...
        'Name', sprintf('Zad.5: %s - przed/po', r.label));

    % Subplot 1: Korelacja krzyżowa
    subplot(2,3,[1 2]);
    ch = 1;
    len_def = out_def.length_domain;
    len_opt = out_opt.length_domain;
    D_s_def = out_def.D_s;
    D_s_opt = out_opt.D_s;

    plot(len_def, out_def.dane_bez_obwiedni(ch,:), 'b-', 'LineWidth', 1.5); hold on;
    plot(len_opt, out_opt.dane_bez_obwiedni(ch,:), 'r-', 'LineWidth', 1.5);
    D_min = min([params_def.fbg.D, params_opt.fbg.D]);
    D_max = max([params_def.fbg.D, params_opt.fbg.D]);
    for j = 1:N_s
        xline(params_def.fbg.D(j), '--k', '', 'LineWidth', 0.8);
    end
    xlim([D_min-30, D_max+30]);
    xlabel('Pozycja [m]'); ylabel('Amplituda');
    title(sprintf('Korelacja krzyżowa \\lambda_1 - %s', r.label));
    legend(sprintf('Default (MAE=%.0f pm)', out_def.MAE), ...
           sprintf('%s opt (MAE=%.0f pm)', best_name, out_opt.MAE), 'Location', 'best');
    grid on; set(gca, 'FontSize', 11);

    % Subplot 2: Błędy per siatka
    subplot(2,3,3);
    bar([out_def.lB_errors'*1000, out_opt.lB_errors'*1000]);
    xlabel('Siatka'); ylabel('\Delta\lambda_B [pm]');
    title('Błędy detekcji per siatka');
    legend('Default', sprintf('%s opt', best_name), 'Location', 'best');
    grid on; set(gca, 'FontSize', 10);
    if N_s <= 10
        xticklabels(arrayfun(@(x) sprintf('%d', x), 1:N_s, 'UniformOutput', false));
    end

    % Subplot 3: Widma - default
    subplot(2,3,4);
    sim_nm = out_def.simulation_range*1e9;
    colors = lines(min(N_s, 5));
    for i = 1:min(N_s, 5)
        plot(sim_nm, out_def.R_simulated{i}, '-', 'Color', colors(i,:), 'LineWidth', 1.2); hold on;
    end
    xlabel('\lambda [nm]'); ylabel('R');
    title(sprintf('Widma DEFAULT (\\Deltan_{eff}=%.1e)', 0.3e-4));
    grid on; set(gca, 'FontSize', 10);

    % Subplot 4: Widma - optimum
    subplot(2,3,5);
    sim_nm_opt = out_opt.simulation_range*1e9;
    for i = 1:min(N_s, 5)
        plot(sim_nm_opt, out_opt.R_simulated{i}, '-', 'Color', colors(i,:), 'LineWidth', 1.2); hold on;
    end
    xlabel('\lambda [nm]'); ylabel('R');
    title(sprintf('Widma %s (\\Deltan_{eff}=%.1e, spread=%.2fnm)', best_name, best.x(1), best.x(5)));
    grid on; set(gca, 'FontSize', 10);

    % Subplot 5: Podsumowanie - bar MAE + parametry
    subplot(2,3,6);
    bar([out_def.MAE, out_opt.MAE], 'FaceColor', 'flat');
    colordata = [0.6 0.6 0.8; 0.2 0.7 0.3];
    b = bar([out_def.MAE, out_opt.MAE], 'FaceColor', 'flat');
    b.CData = colordata;
    set(gca, 'XTickLabel', {'Default', sprintf('%s opt', best_name)});
    ylabel('MAE [pm]');
    title(sprintf('Poprawa: %.0f → %.0f pm (−%.0f%%)', ...
        out_def.MAE, out_opt.MAE, (out_def.MAE-out_opt.MAE)/out_def.MAE*100));
    grid on; set(gca, 'FontSize', 11);

    sgtitle(sprintf('Wariant %d: %s - porównanie przed/po optymalizacji', v, r.label), ...
        'FontSize', 13, 'FontWeight', 'bold');

    savefig(gcf, fullfile(results_dir, sprintf('comparison_v%d.fig', v)));
    saveas(gcf, fullfile(results_dir, sprintf('comparison_v%d.png', v)));
end

%% ===== TABELA ZBIORCZA =====
fprintf('\n\n============================================\n');
fprintf('  TABELA ZBIORCZA: PRZED vs PO OPTYMALIZACJI\n');
fprintf('============================================\n\n');

fprintf('%-25s | %-10s | %-10s | %-10s | %-10s\n', ...
    'Wariant', 'Default', 'GA', 'PSO', 'Poprawa');
fprintf('%s\n', repmat('-', 1, 75));
for v = 1:length(all_results)
    r = all_results(v);
    best_mae = min(r.ga.mae, r.pso.mae);
    imp = (r.default_mae - best_mae) / r.default_mae * 100;
    fprintf('%-25s | %-10.1f | %-10.1f | %-10.1f | −%.0f%%\n', ...
        r.label, r.default_mae, r.ga.mae, r.pso.mae, imp);
end

%% ===== WYKRES ZBIORCZY =====
figure('Position', [50 50 1000 500], 'theme', 'light', 'Name', 'Zad.5: Podsumowanie');

mae_data = zeros(length(all_results), 3);
labels_v = {};
for v = 1:length(all_results)
    mae_data(v,:) = [all_results(v).default_mae, all_results(v).ga.mae, all_results(v).pso.mae];
    labels_v{v} = all_results(v).label;
end

b = bar(mae_data);
b(1).FaceColor = [0.7 0.7 0.7];
b(2).FaceColor = [0.2 0.5 0.8];
b(3).FaceColor = [0.8 0.3 0.2];
set(gca, 'XTickLabel', labels_v, 'XTickLabelRotation', 20);
ylabel('MAE [pm]');
title(sprintf('Porównanie: Default vs GA vs PSO - %s', datestr(now, 'yyyy-mm-dd HH:MM')));
legend('Default (p=8, \Deltan_{eff}=0.3e-4)', 'GA optimum', 'PSO optimum', 'Location', 'northwest');
grid on; set(gca, 'FontSize', 12);

savefig(gcf, fullfile(results_dir, 'comparison_summary.fig'));
saveas(gcf, fullfile(results_dir, 'comparison_summary.png'));

%% ===== ZAPIS =====
save(fullfile(results_dir, 'comparison_results.mat'), 'all_results');
fprintf('\nWszystko zapisane w %s\n', results_dir);
fprintf('Skrypt zakończony: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%% ===== HELPER =====
function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end
