function [lambdas] = generateLambdas(laser)
% GENERATELAMBDAS Generate wavelength array for FBG interrogation
%
% [lambdas] = generateLambdas(laser)
%
% Inputs:
%   laser - Structure with laser parameters:
%       laser.lasers_mode - 'pm' or 'number'
%       laser.step - Step size (in pm) or number of wavelengths
%       laser.wavelength_range - [min_wavelength, max_wavelength] in nm
%
% Outputs:
%   lambdas - Array of wavelengths in meters

% Check laser parameters
if ~contains(laser.lasers_mode, {'pm', 'number'})
    error('Invalid laser mode: must be "pm" or "number"');
end

if numel(laser.wavelength_range) ~= 2 || ~isnumeric(laser.wavelength_range)
    error('Wavelength range must be a numeric array with two elements');
end

% Sort wavelength range
laser.wavelength_range = sort(laser.wavelength_range);

% Check wavelength range validity (typically C-band or L-band)
if laser.wavelength_range(1) < 1200 || laser.wavelength_range(1) > 1600
   error('Wavelength range outside typical telecom bands (1200-1600 nm)'); 
end

% Generate wavelength array based on mode
switch laser.lasers_mode
    case 'pm'
        % Fixed step in picometers
        lambdas = laser.wavelength_range(1):laser.step*1e-3:laser.wavelength_range(2);
    case 'number'
        % Fixed number of wavelengths
        lambdas = linspace(laser.wavelength_range(1), laser.wavelength_range(2), laser.step);
end

% Convert from nm to m
lambdas = lambdas * 1e-9;
end
