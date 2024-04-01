function [simulation_range, spectra] = generateLasers(laser)
% GENERATELASERS Generate laser source spectral profiles
%
% [simulation_range, spectra] = generateLasers(laser)
%
% Inputs:
%   laser - Structure with laser parameters:
%       laser.lasers_mode - 'pm' or 'number'
%       laser.step - Step size (in pm) or number of wavelengths
%       laser.wavelength_range - [min_wavelength, max_wavelength] in nm
%       laser.FWHM - Full width at half maximum (in pm)
%       laser.power - Optical power (in W)
%       laser.shape - 'lorentz' or 'gauss'
%
% Outputs:
%   simulation_range - Wavelength range for simulation (in m)
%   spectra - Spectral distribution of each laser

% Generate wavelength array for laser center wavelengths
switch laser.lasers_mode
    case 'pm'
        lambdas = laser.wavelength_range(1):laser.step*1e-3:laser.wavelength_range(2);
    case 'number'
        lambdas = linspace(laser.wavelength_range(1), laser.wavelength_range(2), laser.step);
end

% Convert from nm to m
lambdas = lambdas * 1e-9;

% Generate simulation range with finer resolution around laser wavelengths
fwhm_m = laser.FWHM * 1e-12;  % Convert pm to m
simulation_range = (lambdas(1) - 2*fwhm_m):(1e-12):(lambdas(end) + 2*fwhm_m);

% Generate laser spectral profiles
spectra = laserModel(simulation_range, lambdas, laser);
end

function spectra = laserModel(lambdas, center_lambda, laser)
% Generate spectral profiles for each laser wavelength

% Convert FWHM from pm to m
FWHM = laser.FWHM * 1e-12;

% Initialize spectra matrix
spectra = zeros(length(center_lambda), length(lambdas));

% Generate spectrum for each center wavelength
for i = 1:length(center_lambda)
    switch laser.shape
        case 'gauss'
            % Gaussian spectral profile
            c = FWHM / (2 * sqrt(2 * log(2)));
            a = laser.power / (c * sqrt(2 * pi));
            spectra(i, :) = a * exp(-(lambdas - center_lambda(i)).^2 ./ (2 * c^2));
            
        case 'lorentz'
            % Lorentzian spectral profile
            spectra(i, :) = laser.power * 1/pi * (0.5 * FWHM) ./ ...
                           ((lambdas - center_lambda(i)).^2 + (0.5 * FWHM)^2);
            
        otherwise
            error('Unknown laser shape! Supported shapes: gauss, lorentz');
    end
end
end
