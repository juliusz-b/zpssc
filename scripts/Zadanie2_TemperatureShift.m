%% Zadanie2_TemperatureShift.m
% Demonstracja wpływu przesunięcia temperaturowego λ_B na detekcję
%
% Projekt: PN/01/0321/2022

clear all; close all; clc;
cd('C:\Users\jbojarczuk\OneDrive - Politechnika Warszawska\PRACA\Insytut Telekomunikacji\Perly Nauki 2022\zpssc');
AddAllSubfolders();

%% ===== WARIANTY =====
shift_variants = [0, 200e-12, 500e-12]; % max shift [m]
variant_labels = {'Baseline (0 pm)', 'Gradient 0→200pm', 'Gradient 0→500pm'};

results = {};
for v = 1:length(shift_variants)
    params = defaultParams();
    params.show_plots = false;
    params.fbg.lambda_shifts = linspace(0, shift_variants(v), params.fbg.N_s);

    fprintf('Symulacja: %s... ', variant_labels{v}); tic;
    results{v} = runSimulation(params);
    fprintf('MAE=%.1f pm, t=%.1fs\n', results{v}.MAE, toc);
end

%% ===== WYNIKI TEKSTOWE =====
fprintf('\n========================================\n');
fprintf('  WYNIKI PRZESUNIĘCIA TEMPERATUROWEGO\n');
fprintf('========================================\n');
for v = 1:length(shift_variants)
    out = results{v};
    fprintf('\n--- %s ---\n', variant_labels{v});
    for i = 1:length(out.lB_org)
        fprintf('  FBG%d: org=%.4f nm, rec=%.4f nm, err=%+.1f pm\n', ...
            i, out.lB_org(i), out.lB_rec(i), out.lB_errors(i)*1000);
    end
    fprintf('  MAE = %.1f pm\n', out.MAE);
end

%% ===== WYKRESY =====

% Widma FBG
figure('Position', [50 50 1400 500], 'theme', 'light', 'Name', 'Zad.2: Widma FBG');
colors = lines(5);
for v = 1:length(shift_variants)
    subplot(1, length(shift_variants), v);
    sim_nm = results{v}.simulation_range * 1e9;
    for i = 1:length(results{v}.R_simulated)
        plot(sim_nm, results{v}.R_simulated{i}, '-', 'Color', colors(i,:), 'LineWidth', 1.5);
        hold on;
    end
    xlabel('\lambda [nm]'); ylabel('Reflektancja');
    title(variant_labels{v});
    grid on; set(gca, 'FontSize', 11);
end

% Błędy detekcji
figure('Position', [50 600 1200 500], 'theme', 'light', 'Name', 'Zad.2: Błędy detekcji');
N_s = length(results{1}.lB_org);

subplot(1,2,1);
err_mat = zeros(N_s, length(shift_variants));
for v = 1:length(shift_variants)
    err_mat(:,v) = results{v}.lB_errors' * 1000;
end
bar(err_mat);
xlabel('Numer siatki'); ylabel('Błąd \Delta\lambda_B [pm]');
title('Błędy detekcji per siatka');
legend(variant_labels, 'Location', 'best');
grid on; set(gca, 'FontSize', 11);
xticklabels(arrayfun(@(x) sprintf('FBG%d', x), 1:N_s, 'UniformOutput', false));

subplot(1,2,2);
mae_vals = cellfun(@(r) r.MAE, results);
bar(mae_vals, 'FaceColor', [0.2 0.6 0.8]);
ylabel('MAE [pm]'); title('MAE per wariant');
set(gca, 'XTickLabel', variant_labels, 'XTickLabelRotation', 20);
grid on; set(gca, 'FontSize', 11);

% Scatter: oryginalne vs odtworzone
figure('Position', [50 50 800 600], 'theme', 'light', 'Name', 'Zad.2: lambda_B scatter');
markers = {'o', 's', 'd'};
for v = 1:length(shift_variants)
    plot(results{v}.lB_org, results{v}.lB_rec, markers{v}, 'MarkerSize', 10, ...
        'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
    hold on;
end
all_org = cellfun(@(r) r.lB_org, results, 'UniformOutput', false);
all_org = [all_org{:}];
lim = [min(all_org)-0.05, max(all_org)+0.05];
plot(lim, lim, 'k--', 'LineWidth', 1);
xlabel('Oryginalna \lambda_B [nm]'); ylabel('Odtworzona \lambda_B [nm]');
title('Dokładność odtwarzania \lambda_B');
legend([variant_labels, {'Idealna'}], 'Location', 'southeast');
grid on; set(gca, 'FontSize', 12);
axis equal; xlim(lim); ylim(lim);

fprintf('\nWykresy gotowe.\n');
