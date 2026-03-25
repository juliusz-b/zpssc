%% WP2_SystemSimulation_v2.m
% Główna symulacja systemu czujnikowego FBG z multipleksacją kodową (CDM)
% Używa runSimulation() - czysta funkcja bez plotów, wyniki w workspace.
%
% Projekt: PN/01/0321/2022

clear all; close all; clc;
cd('C:\Users\jbojarczuk\OneDrive - Politechnika Warszawska\PRACA\Insytut Telekomunikacji\Perly Nauki 2022\zpssc');
AddAllSubfolders();

%% ===== PARAMETRY =====
params = defaultParams();
params.show_plots = true;  % wyświetl wykresy

% Ewentualne modyfikacje:
% params.p = 10;
% params.fbg.max_bounces = 3;
% params.fbg.deltaneffs = 1e-4 * ones(1, 5);
% params.fbg.lambda_shifts = linspace(0, 0.2e-9, 5);

%% ===== SYMULACJA =====
fprintf('Uruchamiam symulację: %s, %s, p=%d, N_s=%d\n', ...
    params.scenariusz, params.kod, params.p, params.fbg.N_s);

tic;
out = runSimulation(params);
fprintf('Czas: %.2f s\n', toc);

%% ===== WYNIKI =====
fprintf('\n========================================\n');
fprintf('  WYNIKI DETEKCJI lambda_B\n');
fprintf('========================================\n');
fprintf('%-10s %-16s %-16s %-12s\n', 'Siatka', 'Oryginalne [nm]', 'Odtworzone [nm]', 'Błąd [pm]');
fprintf('----------------------------------------\n');
for i = 1:params.fbg.N_s
    fprintf('FBG %-6d %-16.4f %-16.4f %-12.2f\n', ...
        i, out.lB_org(i), out.lB_rec(i), out.lB_errors(i)*1000);
end
fprintf('----------------------------------------\n');
fprintf('MAE  = %.2f pm\n', out.MAE);
fprintf('MAX  = %.2f pm\n', out.MAX_err);
