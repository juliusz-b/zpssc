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
%       fbg.model - (optional) 'tanh' (default, fast) or 'tmm' (transfer matrix)
%       fbg.lambda_shifts - (optional) Bragg wavelength shifts [m] per grating
%                           Default: zeros (no shift). Positive = red shift.
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
sections = 100;      % Number of sections for TMM simulation

% Select model
if isfield(fbg, 'model')
    model = fbg.model;
else
    model = 'tanh';
end

% Apply wavelength shifts by modifying periods
% lambda_B = 2 * neff * period => period_shifted = period * (1 + shift/lambda_B)
periods = fbg.periods;
if isfield(fbg, 'lambda_shifts') && any(fbg.lambda_shifts ~= 0)
    has_shifts = true;
    for i = 1:fbg.N_s
        lambda_B_base = 2 * neff * periods(i);
        periods(i) = periods(i) * (1 + fbg.lambda_shifts(i) / lambda_B_base);
    end
else
    has_shifts = false;
end

% Initialize output cell arrays
R_s = cell(1, fbg.N_s);
R_simulated = cell(1, fbg.N_s);

% Generate reflectance spectra for each grating at the given wavelengths
for i = 1:fbg.N_s
    R_s{i} = generateFBG(lambdas, periods(i), neff, ...
        fbg.deltaneffs(i), fbg.grating_lengths(i), visibility, sections, model);
end

% If a detailed simulation range is provided, calculate higher-resolution spectra
if nargin > 2
    % Check if all gratings are identical (optimize only if no shifts)
    all_identical = ~has_shifts && ...
        length(unique(fbg.deltaneffs)) == 1 && ...
        length(unique(periods)) == 1 && ...
        length(unique(fbg.grating_lengths)) == 1;

    if all_identical
        % Calculate only one spectrum and replicate it
        R_simulated{1} = generateFBG(simulation_range, periods(1), ...
            neff, fbg.deltaneffs(1), fbg.grating_lengths(1), visibility, sections, model);
        for i = 2:fbg.N_s
            R_simulated{i} = R_simulated{1};
        end
    else
        % Calculate spectra for each grating individually
        for i = 1:fbg.N_s
            R_simulated{i} = generateFBG(simulation_range, periods(i), ...
                neff, fbg.deltaneffs(i), fbg.grating_lengths(i), visibility, sections, model);
        end
    end
end
end

function R = generateFBG(lambdas, period, neff, deltaneff, grating_length, visibility, sections, model)
% Helper: dispatch to the selected FBG model
    switch model
        case 'tanh'
            R = generateBraggGratingTanh(lambdas, period, neff, deltaneff, grating_length, visibility);
        case 'tmm'
            R = generateBraggGrating(lambdas, period, neff, deltaneff, grating_length, visibility, sections);
        otherwise
            error('Unknown FBG model: %s. Use ''tanh'' or ''tmm''.', model);
    end
end
