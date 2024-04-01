function plotGratingsSnr(lambdas, R_s, SNR_ch, fbg)
% PLOTGRATINGSSNR Visualize SNR for FBG array
%
% plotGratingsSnr(lambdas, R_s, SNR_ch, fbg)
%
% Inputs:
%   lambdas - Wavelength array [m]
%   R_s - Cell array with FBG reflectance spectra
%   SNR_ch - Channel SNR [dB]
%   fbg - (optional) FBG parameters structure

% Check for sufficient data
if length(R_s) < 2
    warning('Insufficient gratings to plot SNR');
    return;
end

% Number of gratings
N_s = length(R_s);

% Convert cell array to matrix
step = length(lambdas);
R_s_mat = cell2mat(R_s);
R = [];
for i = 1:N_s
    R = [R; R_s_mat(1:step)];
    R_s_mat(1:step) = [];
end

% Small epsilon to avoid division by zero
epsilon = 0;

% Initialize SNR matrix
SN = zeros(length(lambdas), N_s);
Rk = [];

% Calculate total reflection
RTT = 0;
for i = 1:length(lambdas)
    RTT = RTT + R(1, i);
    
    % Add reflections from other gratings
    if N_s > 1
        rs = 0;
        for j = 2:N_s
            rr = 1;
            for k = 1:j-1
                rr = rr * (1-R(k, i))^2;
            end
            rs = rs + rr * R(j, i);
        end
        RTT = RTT + rs;
    end
end

% Calculate SNR for each wavelength and grating
for alfa = 1:length(lambdas)
    % First grating
    Rk(1, alfa) = R(1, alfa);
    SN(alfa, 1) = Rk(1, alfa) / (RTT - Rk(1, alfa) + 10^(-SNR_ch/10) + epsilon);
    
    % Other gratings
    if N_s > 1
        for k = 1:N_s-1
            rr = 1;
            for i = 1:k-1
                rr = rr * (1-R(i, alfa))^2;
            end
            
            Rk(k, alfa) = rr * R(k, alfa);
            SN(alfa, k) = Rk(k, alfa) / (RTT - Rk(k, alfa) + 10^(-SNR_ch/10) + epsilon);
        end
    end
end

% Convert to dB scale
SN = real(10 * log10(SN));

% Create figure
figure('Renderer', 'painters', 'Position', [100 100 900 600]);
surf(1:N_s, lambdas*1e9, SN, 'edgecolor', 'none');
alpha(0.7);

% Add labels and colorbar
xlabel('Grating number [#]');
ylabel('Wavelength [nm]');
zlabel('SNR_{receiver} [dB]');

% Add reference lines
hold on;
plot3(repmat([1], length(lambdas), 1), lambdas*1e9, repmat(SNR_ch, length(lambdas), 1), 'linewidth', 3, 'color', 'r');
plot3(repmat([1], length(lambdas), 1), lambdas*1e9, repmat(max(SN, [], 'all'), length(lambdas), 1), 'linewidth', 3, 'color', [0.8 0 0]);

% Add arrows to show SNR range
p1 = [1 lambdas(round(length(lambdas)/2))*1e9 SNR_ch];
p2 = [1 lambdas(round(length(lambdas)/2))*1e9 max(SN, [], 'all')];
line([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);

% Add original FBG spectrum if fbg structure provided
if nargin > 3
    lmb = linspace(lambdas(1), lambdas(end), 1000);
    temp_N_s = fbg.N_s;
    fbg.N_s = 1;
    Rmb = generateAllGratings(lmb, fbg);
    fbg.N_s = temp_N_s;
    plot3(repmat([1], length(lmb), 1), lmb*1e9, 10*log10(Rmb{1}*10^(SNR_ch/10)), 'linewidth', 2);
end

% Set axis limits and colorbar
zlim([min(SN(SN~=-Inf), [], 'all') SNR_ch]);
c = colorbar;
caxis([min(SN, [], 'all') max(SN, [], 'all')]);
c.Label.String = 'SNR_{receiver} [dB]';
colormap(jet(10));
view([130 16]);
set(gca, 'fontsize', 15);
end
