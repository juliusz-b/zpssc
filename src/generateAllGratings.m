function [R_s, R_simulated] = generateAllGratings(lambdas, fbg, simulation_range)
% GENERATEALLGRATINGS Generate multiple FBGs with the given parameters
%
% [R_s, R_simulated] = generateAllGratings(lambdas, fbg, simulation_range)
%
% Inputs:
%   lambdas - Wavelength array for sensor interrogation [m]
%   fbg - Structure containing FBG parameters:
%       fbg.N_s - Number of gratings
%       fbg.periods - Period of each grating [m]
%       fbg.deltaneffs - Index modulation depth of each grating
%       fbg.grating_lengths - Length of each grating [m]
%   simulation_range - Wavelength range for detailed simulation [m]
%
% Outputs:
%   R_s - Cell array with reflectance spectrum for each grating at the
%         wavelengths specified in 'lambdas'
%   R_simulated - Cell array with reflectance spectrum for each grating
%                at the wavelengths specified in 'simulation_range'

% Define standard FBG parameters
neff = 1.447;        % Effective refractive index
visibility = 1;      % Fringe visibility
sections = 100;      % Number of sections for simulation

% Initialize output cell arrays
R_s = cell(1, fbg.N_s);
R_simulated = cell(1, fbg.N_s);

% Calculate theoretical maximum reflectance for debugging
lD = 2 * neff * fbg.periods(1);
lmbmax = (1 + fbg.deltaneffs(1) / neff) * lD;
kappa = pi / lmbmax * visibility * fbg.deltaneffs(1);
tanHH = @(x) (exp(x) - exp(-x)) ./ (exp(x) + exp(-x));
rmax = tanHH(kappa * fbg.grating_lengths(1))^2;

% Generate reflectance spectra for each grating at the given wavelengths
parfor i = 1:fbg.N_s
    R_s{i} = generateBraggGrating(lambdas, fbg.periods(i), neff, ...
        fbg.deltaneffs(i), fbg.grating_lengths(i), visibility, sections);
end

% If a detailed simulation range is provided, calculate higher-resolution spectra
if nargin > 2
    % Check if all gratings are identical to optimize computation
    if length(unique(fbg.deltaneffs)) == 1 && ...
       length(unique(fbg.periods)) == 1 && ...
       length(unique(fbg.grating_lengths)) == 1
        
        % Calculate only one spectrum and replicate it
        R_simulated{1} = generateBraggGrating(simulation_range, fbg.periods(1), ...
            neff, fbg.deltaneffs(1), fbg.grating_lengths(1), visibility, sections);
        
        for i = 2:fbg.N_s
            R_simulated{i} = R_simulated{1};
        end
    else
        % Calculate spectra for each grating individually
        parfor i = 1:fbg.N_s
            R_simulated{i} = generateBraggGrating(simulation_range, fbg.periods(i), ...
                neff, fbg.deltaneffs(i), fbg.grating_lengths(i), visibility, sections);
        end
    end
end
end
