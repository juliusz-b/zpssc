function plotXcovs(data_xcovs, x_s)
% PLOTXCOVS Plot cross-covariance results
%
% plotXcovs(data_xcovs, x_s)
%
% Inputs:
%   data_xcovs - Cross-covariance data matrix
%   x_s - (optional) X-axis values

% Set default x-axis if not provided
if nargin < 2
    x_s = 0:length(data_xcovs(1, 1:end))-1;
    xlab = 'Delay [samples]';
else
    xlab = 'Delay [m]';
end

% Create legend entries
leg = {};
for i = 1:size(data_xcovs, 1)
    leg{i} = ['Laser ' sprintf('%i', i)];
end

% Create figure with specified size
figure('Renderer', 'painters', 'Position', [100 100 900 600]);

% Plot each cross-covariance
for i = 1:size(data_xcovs, 1)
    plot(x_s, data_xcovs(i, 1:end), 'linewidth', 3);
    hold on;
end

% Set plot labels and properties
xlabel(xlab);
ylabel('XCov amplitude [a.u.]')
ylim([0 (max(data_xcovs(:)))]);
grid minor;
set(gca, 'fontsize', 15);

% Add legend if multiple signals
if length(leg) > 1
    legend(leg, 'Box', 'off', 'Location', 'best');
end
end
