function C = figStyle(fig, widthIn, heightIn)
% FIGSTYLE One visual language for the figures shown in the README.
%
%   C = figStyle()                 returns the palette only
%   C = figStyle(fig, w, h)        sizes the figure to w x h inches and
%                                  restyles every axes and legend in it
%
% The look mirrors python_selfcal_cdm/figstyle.py: sans-serif type, boxed
% axes with major and minor ticks pointing inwards, bold axis labels, framed
% legends, Paul Tol's colour-blind safe "vibrant" palette. Panel letters are
% bold lowercase, added with figPanelLetter().
%
% Palette roles (same quantity, same colour in every figure):
%   C.blue    direct return, CDM, the better configuration
%   C.verm    measured or uncorrected, TDM, the worse configuration
%   C.green   corrected or improved result
%   C.orange  noise, ghosts
%   C.grey    default or reference value
%   C.ltgrey  true (isolated) spectrum, drawn thick
%   C.black   analytical rule or bound

C.blue   = [0 119 187] / 255;
C.sky    = [51 187 238] / 255;
C.orange = [238 119 51] / 255;
C.verm   = [204 51 17] / 255;
C.green  = [0 153 136] / 255;
C.purple = [238 51 119] / 255;
C.grey   = [77 77 77] / 255;
C.ltgrey = [154 154 154] / 255;
C.black  = [34 34 34] / 255;
C.font   = 'Arial';
C.base   = 8;            % base font size in points

if nargin == 0
    return;
end

fig.Units = 'inches';
fig.Position(3:4) = [widthIn, heightIn];
fig.Color = 'w';
fig.PaperPositionMode = 'auto';

ax = findobj(fig, 'Type', 'axes');
for k = 1:numel(ax)
    a = ax(k);
    set(a, 'FontName', C.font, 'FontSize', C.base, ...
        'Box', 'on', 'LineWidth', 0.8, 'TickDir', 'in', ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'TickLength', [0.015 0.01], 'XGrid', 'off', 'YGrid', 'off', ...
        'Layer', 'top', 'XColor', 'k', 'YColor', 'k');
    a.XLabel.FontWeight = 'bold';
    a.YLabel.FontWeight = 'bold';
    a.XLabel.FontSize = C.base + 1;
    a.YLabel.FontSize = C.base + 1;
    a.Title.FontWeight = 'normal';
    a.Title.FontSize = C.base;
    if ~isempty(a.Legend)
        set(a.Legend, 'Box', 'on', 'EdgeColor', 'k', 'Color', 'w', ...
            'FontSize', C.base - 1, 'LineWidth', 0.6);
    end
end

end
