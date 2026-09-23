%% README_Figures.m
% Redraws the simulator figures shown in README.md in one common style
% (src/plots/figStyle.m) from the saved results in results/. No simulation is
% run here, so the script takes seconds. Output: docs/figures/*.png.
%
% Sources:
%   results/WP2/benchmark/noise_vs_sidelobes.fig
%   results/WP2/benchmark/spectrum_reconstruction_p8_vs_p10.fig -> fig_simulator_output.png
%   results/WP3/sensitivity/sensitivity_results.mat             -> fig_sensitivity.png
%   results/WP3/optimization/full_optimization_results.mat
%   results/WP3/optimization/gauss_optimization_results.mat     -> fig_optimization.png

clear; close all;
repoRoot = fileparts(fileparts(mfilename('fullpath')));
cd(repoRoot);
AddAllSubfolders;

outDir = fullfile('docs', 'figures');
if ~exist(outDir, 'dir'), mkdir(outDir); end
C = figStyle();
W = 7.1;   % full-width figure [in]

%% 1. What the simulator returns: correlation trace and a reconstructed spectrum
% (a) correlation vs position with the five gratings marked, NEP 15 pW/sqrt(Hz)
h = openfig(fullfile('results', 'WP2', 'benchmark', 'noise_vs_sidelobes.fig'), 'invisible');
axAll = findobj(h, 'Type', 'axes');
axMain = axAll(2);
withNoise = findobj(axMain, 'Type', 'line', '-not', 'DisplayName', 'Brak szumu');
corr.x = withNoise.XData; corr.y = withNoise.YData;
close(h);

% (b) spectrum of one grating reconstructed with p = 8 and p = 10, plus the true one
h = openfig(fullfile('results', 'WP2', 'benchmark', 'spectrum_reconstruction_p8_vs_p10.fig'), 'invisible');
axAll = findobj(h, 'Type', 'axes');
spec = struct('p', {}, 'k', {}, 'lam', {}, 'rec', {}, 'true', {});
for a = 1:numel(axAll)
    tok = regexp(axAll(a).Title.String, 'p=(\d+), FBG(\d+)', 'tokens', 'once');
    L = findobj(axAll(a), 'Type', 'line');
    s.p = str2double(tok{1}); s.k = str2double(tok{2});
    for l = 1:numel(L)
        if isequal(L(l).Color, [0 0 0]), s.true = L(l).YData; else, s.rec = L(l).YData; end
        s.lam = L(l).XData;
    end
    spec(end+1) = s; %#ok<SAGROW>
end
close(h);
kShow = 3;
s8 = spec([spec.p] == 8 & [spec.k] == kShow);
s10 = spec([spec.p] == 10 & [spec.k] == kShow);

fbgPos = 400:20:480;
scale = 1e4;   % correlation amplitude shown in units of 1e-4

fig = figure('Visible', 'off');
t = tiledlayout(fig, 1, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(t, [1 3]);
plot(ax1, corr.x, corr.y * scale, 'Color', C.blue, 'LineWidth', 0.8); hold(ax1, 'on');
yl = [-3 3];
for k = 1:numel(fbgPos)
    plot(ax1, fbgPos(k) * [1 1], yl, ':', 'Color', C.grey, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    lbl = sprintf('FBG%d', k);
    if k == kShow, lbl = [lbl ' (b)']; end %#ok<AGROW>
    text(ax1, fbgPos(k), yl(2) * 0.95, lbl, 'FontName', C.font, ...
        'FontSize', C.base - 1, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 1);
end
xlim(ax1, [370 510]); ylim(ax1, yl);
xlabel(ax1, 'Position along the fiber [m]');
ylabel(ax1, 'Correlation [\times10^{-4}]');
figPanelLetter(ax1, 'a');

ax2 = nexttile(t, [1 2]);
plot(ax2, s8.lam, s8.true, 'Color', C.ltgrey, 'LineWidth', 2.4); hold(ax2, 'on');
plot(ax2, s8.lam, s8.rec, '-o', 'Color', C.verm, 'LineWidth', 1.0, 'MarkerSize', 3, 'MarkerFaceColor', C.verm);
plot(ax2, s10.lam, s10.rec, '-o', 'Color', C.blue, 'LineWidth', 1.0, 'MarkerSize', 3, 'MarkerFaceColor', C.blue);
xlim(ax2, [1534.3 1536.7]); ylim(ax2, [0 1.08]);
xlabel(ax2, 'Wavelength [nm]');
ylabel(ax2, 'Normalized spectrum');
legend(ax2, {'True spectrum', 'p = 8 (255 chips)', 'p = 10 (1023 chips)'}, 'Location', 'northwest');
figPanelLetter(ax2, 'b');
figStyle(fig, W, 2.5);
exportgraphics(fig, fullfile(outDir, 'fig_simulator_output.png'), 'Resolution', 300);
close(fig);

%% 2. Sensitivity: which parameter matters
S = load(fullfile('results', 'WP3', 'sensitivity', 'sensitivity_results.mat'));
fig = figure('Visible', 'off');
t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(t);
plot(ax, S.p_v, S.mae_p, '-o', 'Color', C.blue, 'MarkerFaceColor', C.blue, 'LineWidth', 1.2, 'MarkerSize', 4);
xlabel(ax, 'Code exponent p (2^p - 1 chips)'); ylabel(ax, 'MAE [pm]');
xlim(ax, [3.5 10.5]); ylim(ax, [0 560]); ax.XTick = 4:2:10;
text(ax, S.p_v, S.mae_p + 28, compose('%.0f', S.mae_p), 'FontName', C.font, ...
    'FontSize', C.base - 1, 'HorizontalAlignment', 'center');
figPanelLetter(ax, 'a');

ax = nexttile(t);
names = {'Code length p', 'Temperature gradient', '\Deltan_{eff}', 'Number of gratings N_s', 'Detector NEP'};
range = [max(S.mae_p) - min(S.mae_p), max(S.mae_gr) - min(S.mae_gr), ...
    max(S.mae_dn) - min(S.mae_dn), max(S.mae_ns) - min(S.mae_ns), max(S.mae_nep) - min(S.mae_nep)];
barh(ax, flip(range), 0.6, 'FaceColor', C.blue, 'EdgeColor', 'none');
ax.YTick = 1:5; ax.YTickLabel = flip(names); ax.YMinorTick = 'off';
xlabel(ax, 'MAE range over the sweep [pm]'); xlim(ax, [0 480]);
text(ax, flip(range) + 10, 1:5, compose('%.0f', flip(range)), 'FontName', C.font, ...
    'FontSize', C.base - 1, 'VerticalAlignment', 'middle');
figPanelLetter(ax, 'b');
figStyle(fig, W, 2.4);
exportgraphics(fig, fullfile(outDir, 'fig_sensitivity.png'), 'Resolution', 300);
close(fig);

%% 3. Optimization: default vs optimized (centroid) vs optimized (Gaussian fit)
F = load(fullfile('results', 'WP3', 'optimization', 'full_optimization_results.mat'));
G = load(fullfile('results', 'WP3', 'optimization', 'gauss_optimization_results.mat'));
nsList = [G.all_gauss.N_s];
vals = zeros(numel(nsList), 3);
for i = 1:numel(nsList)
    r = F.all_results([F.all_results.N_s] == nsList(i) & [F.all_results.grad_nm] == 0);
    vals(i, 1) = r.default_mae;
    vals(i, 2) = min(r.ga.mae, r.pso.mae);
    vals(i, 3) = G.all_gauss(i).pso.mae;
end
fig = figure('Visible', 'off');
ax = axes(fig);
b = bar(ax, vals, 0.75, 'EdgeColor', 'none');
b(1).FaceColor = C.ltgrey; b(2).FaceColor = C.blue; b(3).FaceColor = C.green;
ax.XTick = 1:numel(nsList); ax.XTickLabel = compose('N_s = %d', nsList); ax.XMinorTick = 'off';
ylabel(ax, 'MAE [pm]'); ylim(ax, [0 460]);
for i = 1:numel(nsList)
    for j = 1:3
        text(ax, b(j).XEndPoints(i), vals(i, j) + 10, sprintf('%.1f', vals(i, j)), ...
            'FontName', C.font, 'FontSize', C.base - 1, 'HorizontalAlignment', 'center');
    end
end
legend(ax, {'Default (p = 8, centroid)', 'Optimized, centroid', 'Optimized, Gaussian fit'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal');
figStyle(fig, 4.6, 2.9);
exportgraphics(fig, fullfile(outDir, 'fig_optimization.png'), 'Resolution', 300);
close(fig);

fprintf('README figures written to %s\n', outDir);
