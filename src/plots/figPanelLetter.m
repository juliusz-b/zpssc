function figPanelLetter(ax, letter)
% FIGPANELLETTER Bold lowercase panel letter above the top-left corner.
%
%   figPanelLetter(ax, 'a')
%
% Same convention as figstyle.letter() in python_selfcal_cdm.

text(ax, 0.0, 1.04, letter, 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

end
