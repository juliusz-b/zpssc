%% Zadanie6_BenchmarkCodesDenoising.m
% Benchmark: porównanie kodów (Kasami, Gold) x metod odszumiania x NEP
%
% Projekt: PN/01/0321/2022

clear all; close all; clc;
cd('C:\Users\jbojarczuk\OneDrive - Politechnika Warszawska\PRACA\Insytut Telekomunikacji\Perly Nauki 2022\zpssc');
AddAllSubfolders();

%% ===== KONFIGURACJA BENCHMARKU =====
codes = {'kasami', 'gold'};
denoise_methods = {'none', 'sgolay', 'wavelet'};
NEP_values = [5, 15, 30] * 1e-12;

results_dir = fullfile('results', 'WP2', 'benchmark');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end

%% ===== URUCHOM BENCHMARK =====
n_total = length(codes) * length(denoise_methods) * length(NEP_values);
bench = struct('kod',{}, 'denoise',{}, 'NEP_pW',{}, 'MAE',{}, 'MAX',{}, 'time',{});
idx = 0;

fprintf('=== BENCHMARK: %d konfiguracji ===\n', n_total);
total_tic = tic;

for c = 1:length(codes)
    for d = 1:length(denoise_methods)
        for n = 1:length(NEP_values)
            idx = idx + 1;
            params = defaultParams();
            params.show_plots = false;
            params.kod = codes{c};
            params.denoise_m.method = denoise_methods{d};
            params.noise.NEP = NEP_values(n);

            fprintf('[%2d/%2d] %s+%s+NEP=%.0fpW... ', idx, n_total, ...
                codes{c}, denoise_methods{d}, NEP_values(n)*1e12);
            tic;
            try
                out = runSimulation(params);
                t = toc;
                bench(idx).kod = codes{c};
                bench(idx).denoise = denoise_methods{d};
                bench(idx).NEP_pW = NEP_values(n)*1e12;
                bench(idx).MAE = out.MAE;
                bench(idx).MAX = out.MAX_err;
                bench(idx).time = t;
                fprintf('MAE=%.1f pm, t=%.1fs\n', out.MAE, t);
            catch ME
                toc;
                bench(idx).kod = codes{c};
                bench(idx).denoise = denoise_methods{d};
                bench(idx).NEP_pW = NEP_values(n)*1e12;
                bench(idx).MAE = NaN; bench(idx).MAX = NaN; bench(idx).time = NaN;
                fprintf('BŁĄD: %s\n', ME.message);
            end
        end
    end
end

fprintf('\nCałkowity czas: %.1f s\n', toc(total_tic));

%% ===== WYNIKI =====
fprintf('\n==============================================\n');
fprintf('  BENCHMARK: KODY x ODSZUMIANIE x NEP\n');
fprintf('==============================================\n');
fprintf('%-10s %-10s %-8s %-10s %-10s %-8s\n', 'Kod', 'Denoise', 'NEP[pW]', 'MAE[pm]', 'MAX[pm]', 'Czas[s]');
fprintf('----------------------------------------------\n');
for i = 1:length(bench)
    fprintf('%-10s %-10s %-8.0f %-10.1f %-10.1f %-8.1f\n', ...
        bench(i).kod, bench(i).denoise, bench(i).NEP_pW, bench(i).MAE, bench(i).MAX, bench(i).time);
end

%% ===== ZAPIS =====
save(fullfile(results_dir, 'benchmark_results.mat'), 'bench');
fprintf('\nWyniki zapisane w %s\n', results_dir);

%% ===== WYKRESY =====
figure('Position', [50 50 1000 500], 'theme', 'light', 'Name', 'Zad.6: Benchmark');

% MAE per konfiguracja (NEP=15)
subplot(1,2,1);
idx15 = find([bench.NEP_pW] == 15);
labels15 = arrayfun(@(i) sprintf('%s+%s', bench(i).kod, bench(i).denoise), idx15, 'UniformOutput', false);
mae15 = [bench(idx15).MAE];
bar_colors = repmat([0.2 0.6 0.9], length(idx15), 1);
bar_colors(4:end, :) = repmat([0.9 0.4 0.2], length(idx15)-3, 1);
b = bar(mae15, 'FaceColor', 'flat');
b.CData = bar_colors;
set(gca, 'XTickLabel', labels15, 'XTickLabelRotation', 30);
ylabel('MAE [pm]'); title('MAE przy NEP=15 pW/\surdHz');
grid on; set(gca, 'FontSize', 11);

% MAE vs NEP
subplot(1,2,2);
neps = [5 15 30];
for c = 1:length(codes)
    mae_c = zeros(1, length(neps));
    for n = 1:length(neps)
        idx_cn = find([bench.NEP_pW]==neps(n) & strcmp({bench.kod},codes{c}) & strcmp({bench.denoise},'none'));
        if ~isempty(idx_cn); mae_c(n) = bench(idx_cn(1)).MAE; end
    end
    plot(neps, mae_c, 'o-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'auto');
    hold on;
end
xlabel('NEP [pW/\surdHz]'); ylabel('MAE [pm]');
title('MAE vs NEP (denoise=none)');
legend(codes, 'Location', 'best');
grid on; set(gca, 'FontSize', 11);

savefig(gcf, fullfile(results_dir, 'benchmark_results.fig'));
saveas(gcf, fullfile(results_dir, 'benchmark_results.png'));
fprintf('Wykresy gotowe.\n');
