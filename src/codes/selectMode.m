function [mode_fun, power] = selectMode(mode)
% SELECTMODE Select signal mode for the simulation
%
% [mode_fun, power] = selectMode(mode)
%
% Inputs:
%   mode - 'unipolar' or 'bipolar'
%
% Outputs:
%   mode_fun - Function handle for converting signals
%   power - Power parameter for signal processing

% Convert mode to lowercase for case-insensitive comparison
mode = lower(mode);

% Select appropriate mode function and power parameter
if contains(mode, 'unipolar')
    mode_fun = @toUnipolar;
    power = 1;
else
    mode_fun = @toBipolar;
    power = 1;
end
end

function [x_out] = toUnipolar(x)
% Convert bipolar signal to unipolar
x_out = x;
x_out(x == -1) = 0;
end

function [x_out] = toBipolar(x)
% Convert unipolar signal to bipolar
x_out = x;
x_out(x == 0) = -1;
end
