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
codes.ooc = genSpreadCodes(13, U, 'ooc',[4,2]);
codes.sidelnikov = genSpreadCodes(257, U, 'sidelnikov',[2]);
codes.golay = genSpreadCodes(L, U, 'golay');
codes.chaotic = genSpreadCodes(L, U, 'chaotic');

%% Compare auto-correlation properties
lineColors = [
    0, 0, 1;           % Blue
    1, 0, 0;           % Red
    0, 0.7, 0;         % Green
    0.7, 0, 0.7;       % Purple
    1, 0.5, 0;         % Orange
    0, 0.7, 0.7;       % Cyan
    0.6, 0.3, 0;       % Brown
    0, 0, 0;           % Black
    0.5, 0.5, 0.5;     % Gray
    0.8, 0.6, 0.2      % Gold
    ];
figure('Name', 'Auto-correlation Comparison', 'Position', [100, 100, 1000, 600]);

for i = 1:2
    subplot(1,2, i);

    % Calculate and plot auto-correlation for each code type
    hold on;
    code_types = fieldnames(codes);
    for j = 1:length(code_types)
        code_type = code_types{j};
        code = codes.(code_type)(i, :);

        % Calculate auto-correlation
        [ac,lags] = xcorr(code, code);

        % Plot
        plot(lags,ac / max(ac),'o-' , 'LineWidth', 1.5,'color',lineColors(j,:));
    end

    title(['User ' num2str(i) ' Auto-correlation']);
    xlabel('Lag');
    ylabel('Correlation');

    if i == 1
        legend(code_types, 'Location', 'best','NumColumns',2,'Box','off');
    end

    xlim([-3 3])
    ylim([0 1])

    grid on;
    hold off;
end

%% Compare cross-correlation properties
figure('Name', 'Cross-correlation Matrix', 'Position', [100, 100, 1200, 800],'Color','w');

% Calculate cross-correlation metrics for each code family
cc_metrics = struct();

for i = 1:length(code_types)
    code_type = code_types{i};
    code_set = codes.(code_type);

    % Create cross-correlation matrix
    cc_matrix = ones(U, U);
    max_cc = ones(U, U);

    cross_psnrs = zeros(U,U);
    self_psnrs = zeros(1,U);
    for u1 = 1:U
        self_psnrs(u1) = myPsnr(xcov(code_set(u1, :)));
        for u2 = 1:U
            if u1 ~= u2
                cc = xcov(code_set(u1, :), code_set(u2, :)) / length(code_set(u1, :));
                max_cc(u1, u2) = max(abs(cc));
                cc_matrix(u1, u2) = max(abs(cc));
                cross_psnrs(u1) = myPsnr(cc);
            end
        end
    end

    % Store metrics
    cc_metrics.(code_type).matrix = cc_matrix;
    cc_metrics.(code_type).max = max(max_cc(:));
    cc_metrics.(code_type).mean = mean(max_cc(:));
    cc_metrics.(code_type).std = std(max_cc(:));
    cc_metrics.(code_type).self_psnrs = self_psnrs;
    cc_metrics.(code_type).cross_psnrs = cross_psnrs;

    % Plot the cross-correlation matrix
    subplot(3, 3, i);
    %imagesc(log10(cc_matrix));
    imagesc((cc_matrix));
    colormap(gray);
    colorbar;
    clim([0 1])

    title(['Korelacja wzajemna ' upper(code_type(1:end)) '']);
    xlabel('\lambda_#');
    ylabel('\lambda_#');
    axis square;

    %Add text to display metrics
    %text(0.5, -0.1, ['Max: ' num2str(cc_metrics.(code_type).max, '%.4f')], ...
    %    'Units', 'normalized', 'HorizontalAlignment', 'center');
    % text(0.5, -0
end

%%
lineColors = [
    0, 0, 1;           % Blue
    1, 0, 0;           % Red
    0, 0.7, 0;         % Green
    0.7, 0, 0.7;       % Purple
    1, 0.5, 0;         % Orange
    0, 0.7, 0.7;       % Cyan
    0.6, 0.3, 0;       % Brown
    0, 0, 0;           % Black
    0.5, 0.5, 0.5;     % Gray
    0.8, 0.6, 0.2      % Gold
    ];
figure('color','w');
for i = 1:length(code_types)
    code_type = code_types{i};
    subplot(2,1,1)
    plot(1:U,cc_metrics.(code_type).self_psnrs,'-o','Color',lineColors(i,:));
    hold on;
    xlabel('\lambda_#');
    ylabel('PSNR [dB]');
    title('Auto PSNR')
    subplot(2,1,2)
    plot(1:U,cc_metrics.(code_type).cross_psnrs,'-o','Color',lineColors(i,:));
    hold on;
    xlabel('\lambda_#');
    ylabel('PSNR [dB]');
    title('Cross PSNR')
end
subplot(2,1,1)
legend(upper(code_types),'Box','off','NumColumns',2,'Location','best')

%% filtrowanie

fltr = @(syg,a)  lowpass(syg,a);
a_s = 0.1:0.1:1;

autocors = {};
xcors = {};
figure;
for i = 1:length(code_types)
    code_type = code_types{i};

    for j=1:length(a_s)
        a= a_s(j);
        if a>=1
            syg1 = codes.(code_type)(1,:);
            syg2 = codes.(code_type)(2,:);
        else
            syg1 = fltr(codes.(code_type)(1,:),a);
            syg2 = fltr(codes.(code_type)(2,:),a);
        end
        autocors{i}(j,:) = xcov(syg1);
        xcors{i}(j,:) = xcov(syg1,syg2);
    end
end
%%
lineColors = [
    0, 0, 1;           % Blue
    1, 0, 0;           % Red
    0, 0.7, 0;         % Green
    0.7, 0, 0.7;       % Purple
    1, 0.5, 0;         % Orange
    0, 0.7, 0.7;       % Cyan
    0.6, 0.3, 0;       % Brown
    0, 0, 0;           % Black
    0.5, 0.5, 0.5;     % Gray
    0.8, 0.6, 0.2      % Gold
    ];
figure('color','w');

for i = 1:length(code_types)
    subplot(2,1,1)
    plot(a_s,myPsnr(autocors{i}')-myPsnr(autocors{i}(end,:)),'-o','Color',lineColors(i,:))
    hold on;
xlabel('2f_{pass}/f_s')
    ylabel('\Delta PSNR [dB]');
    title('Auto PSNR')
    subplot(2,1,2)
    plot(a_s,myPsnr(xcors{i}')-myPsnr(xcors{i}(end,:)),'-o','Color',lineColors(i,:))
    hold on;
    %ylim([5 20])
xlabel('2f_{pass}/f_s')
    ylabel('\Delta PSNR [dB]');
    title('Cross PSNR')
end
subplot(2,1,1)
legend(upper(code_types),'Box','off','NumColumns',2,'Location','best')

