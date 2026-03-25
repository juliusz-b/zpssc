%% Zadanie_TDMvsCDM_Ghosty.m
% Porównanie wpływu ghostów wielokrotnych odbić w TDM vs CDM
%
% Cel: wykazać że CDM zapewnia naturalną ochronę przed ghostingiem
%
% Projekt: PN/01/0321/2022
% Data: 2025-12-15

clear all; close all; clc;
cd('C:\Users\jbojarczuk\OneDrive - Politechnika Warszawska\PRACA\Insytut Telekomunikacji\Perly Nauki 2022\zpssc');
AddAllSubfolders();
SILENT_MODE = true;
neff = 1.447;

results_dir = fullfile('results', 'WP2', 'benchmark');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end

%% ===== PARAMETRY =====
dn_vals = [0.3e-4, 0.5e-4, 0.7e-4, 1e-4]; % różne R
scenarios = {'kody2', 'tof'}; % CDM, TDM
bounces_vals = [1, 3];

%% ===== SYMULACJE =====
fprintf('============================================\n');
fprintf('  TDM vs CDM: WPŁYW GHOSTÓW\n');
fprintf('  %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf('============================================\n\n');

mae_results = zeros(length(dn_vals), length(scenarios), length(bounces_vals));

for idn = 1:length(dn_vals)
    dn = dn_vals(idn);
    kappa = pi*dn/(2*neff*530.73e-9);
    Rpk = tanh(kappa*3e-3)^2*100;

    for isc = 1:length(scenarios)
        for ib = 1:length(bounces_vals)
            params = defaultParams();
            params.show_plots = false;
            params.scenariusz = scenarios{isc};
            params.fbg.deltaneffs = dn*ones(1,params.fbg.N_s);
            params.fbg.max_bounces = bounces_vals(ib);

            out = runSimulation(params);
            mae_results(idn, isc, ib) = out.MAE;
        end

        label = upper(strrep(scenarios{isc}, 'kody2', 'cdm'));
        label = strrep(label, 'TOF', 'TDM');
        fprintf('R=%.1f%% %s: 1b=%.1f, 3b=%.1f, ΔMAE=%+.1f pm\n', ...
            Rpk, label, mae_results(idn,isc,1), mae_results(idn,isc,2), ...
            mae_results(idn,isc,2)-mae_results(idn,isc,1));
    end
    fprintf('\n');
end

%% ===== WYKRESY =====
R_vals = zeros(size(dn_vals));
for i=1:length(dn_vals)
    kp=pi*dn_vals(i)/(2*neff*530.73e-9);
    R_vals(i)=tanh(kp*3e-3)^2*100;
end

figure('Position', [50 50 1200 500], 'theme', 'light', 'Name', 'TDM vs CDM: ghosty');

% Panel 1: MAE vs R dla obu systemów
subplot(1,2,1);
plot(R_vals, mae_results(:,1,1), 'bo-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'b'); hold on;
plot(R_vals, mae_results(:,1,2), 'b^--', 'LineWidth', 2, 'MarkerSize', 10);
plot(R_vals, mae_results(:,2,1), 'rs-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(R_vals, mae_results(:,2,2), 'rv--', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('R_{peak} [%]'); ylabel('MAE [pm]');
title('MAE vs reflektancja: TDM vs CDM');
legend('CDM 1-bounce', 'CDM 3-bounce', 'TDM 1-bounce', 'TDM 3-bounce', 'Location', 'best');
grid on; set(gca, 'FontSize', 12);

% Panel 2: Degradacja od ghostów
subplot(1,2,2);
delta_cdm = mae_results(:,1,2) - mae_results(:,1,1);
delta_tdm = mae_results(:,2,2) - mae_results(:,2,1);
bar([delta_cdm, delta_tdm]);
set(gca, 'XTickLabel', arrayfun(@(r) sprintf('R=%.0f%%', r), R_vals, 'UniformOutput', false));
ylabel('\DeltaMAE [pm] (3b - 1b)');
title('Degradacja od ghostów');
legend('CDM', 'TDM', 'Location', 'best');
grid on; set(gca, 'FontSize', 12);

savefig(gcf, fullfile(results_dir, 'tdm_vs_cdm_ghosty.fig'));
saveas(gcf, fullfile(results_dir, 'tdm_vs_cdm_ghosty.png'));

%% ===== PODSUMOWANIE =====
fprintf('============================================\n');
fprintf('  PODSUMOWANIE TDM vs CDM\n');
fprintf('============================================\n\n');
fprintf('%-8s %-8s %-10s %-10s %-12s\n', 'R[%%]', 'System', 'MAE 1b', 'MAE 3b', 'ΔMAE ghosty');
for idn = 1:length(dn_vals)
    for isc = 1:2
        label = {'CDM','TDM'};
        fprintf('%-8.1f %-8s %-10.1f %-10.1f %-+12.1f\n', R_vals(idn), label{isc}, ...
            mae_results(idn,isc,1), mae_results(idn,isc,2), ...
            mae_results(idn,isc,2)-mae_results(idn,isc,1));
    end
end

save(fullfile(results_dir, 'tdm_vs_cdm_results.mat'), ...
    'mae_results', 'dn_vals', 'R_vals', 'scenarios');
fprintf('\nZapisano w %s\n', results_dir);
