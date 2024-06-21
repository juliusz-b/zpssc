%% WP1_CodeAnalysis.m
% This script analyzes various code families for their correlation properties
% and compares their suitability for FBG sensor interrogation systems.

% Initialize workspace
clear all;
close all;
clc;

% Ensure all folders are in the path
AddAllSubfolders;

%% Parameters
p = 8;          % Power of 2 for code length
L = 2^p-1;      % Code length
Nb = 4;         % Samples per bit
U = 8;          % Number of users/codes

%% Generate different code families
codes = struct();
codes.gold = genSpreadCodes(L, U, 'gold');
codes.kasami = genSpreadCodes(L, U, 'kasami');
codes.prbs = genSpreadCodes(L, U, 'prbs');
codes.walsh = genSpreadCodes(L, U, 'walsh');
codes.randi = genSpreadCodes(L, U, 'randi');
codes.ooc = genSpreadCodes(13, 2, 'ooc',[4,2]);
codes.sidelnikov = genSpreadCodes(257, U, 'sidelnikov',[2]);
codes.golay = genSpreadCodes(L, U, 'golay');
codes.chaotic = genSpreadCodes(L, U, 'chaotic');

%% Compare auto-correlation properties
figure('Name', 'Auto-correlation Comparison', 'Position', [100, 100, 1000, 600]);

for i = 1:U
    subplot(2, 4, i);

    % Calculate and plot auto-correlation for each code type
    hold on;
    code_types = fieldnames(codes);
    for j = 1:length(code_types)
        code_type = code_types{j};
        code = codes.(code_type)(i, :);

        % Calculate auto-correlation
        [ac,lags] = xcorr(code, code);

        % Plot
        plot(lags,ac / max(ac) , 'LineWidth', 1.5);
    end

    title(['User ' num2str(i) ' Auto-correlation']);
    xlabel('Lag');
    ylabel('Correlation');

    if i == 1
        legend(code_types, 'Location', 'northeast');
    end

    xlim([-2 2])
    ylim([0 1])

    grid on;
    hold off;
end

%% Compare cross-correlation properties
figure('Name', 'Cross-correlation Matrix', 'Position', [100, 100, 1200, 800]);

% Calculate cross-correlation metrics for each code family
cc_metrics = struct();

for i = 1:length(code_types)
    code_type = code_types{i};
    code_set = codes.(code_type);

    % Create cross-correlation matrix
    cc_matrix = ones(U, U);
    max_cc = ones(U, U);

    for u1 = 1:U
        for u2 = 1:U
            if u1 ~= u2
                cc = xcorr(code_set(u1, :), code_set(u2, :)) / length(code_set(u1, :));
                max_cc(u1, u2) = max(abs(cc));
                cc_matrix(u1, u2) = max(abs(cc));
            end
        end
    end

    % Store metrics
    cc_metrics.(code_type).matrix = cc_matrix;
    cc_metrics.(code_type).max = max(max_cc(:));
    cc_metrics.(code_type).mean = mean(max_cc(:));
    cc_metrics.(code_type).std = std(max_cc(:));

    % Plot the cross-correlation matrix
    subplot(2, 3, i);
    %imagesc(log10(cc_matrix));
    imagesc((cc_matrix));
    colormap(gray);
    colorbar;
    clim([0 1])
    
    title([upper(code_type(1)) code_type(2:end) ' Cross-correlation']);
    xlabel('User Index');
    ylabel('User Index');
    axis square;

    % Add text to display metrics
    % text(0.5, -0.1, ['Max: ' num2str(cc_metrics.(code_type).max, '%.4f')], ...
    %     'Units', 'normalized', 'HorizontalAlignment', 'center');
    % text(0.5, -0
end