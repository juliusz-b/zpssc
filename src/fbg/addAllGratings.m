function [data_out] = addAllGratings(data_in, N_s, D_s, R_s, L_max, R_simulated, lasers_spectra, simulation_range, fbg, fiber)
% ADDALLGRATINGS Add FBG array reflections to input data
%
% [data_out] = addAllGratings(data_in, N_s, D_s, R_s, L_max, R_simulated, 
%                            lasers_spectra, simulation_range, fbg, fiber)
%
% Inputs:
%   data_in - Input optical signal data
%   N_s - Number of gratings
%   D_s - Grating positions (in samples)
%   R_s - Cell array with reflectance spectrum of each grating
%   L_max - Maximum output length
%   R_simulated - Cell array with high-resolution reflectance spectra
%   lasers_spectra - Laser spectral distribution
%   simulation_range - Wavelength range for simulation
%   fbg - FBG parameters structure
%   fiber - Fiber parameters structure
%
% Outputs:
%   data_out - Output optical signal with FBG reflections

% Check input parameters
if length(R_s) < N_s
    error('Insufficient reflectance spectra provided');
end

if length(D_s) < N_s
    error('Insufficient grating positions provided');
end

% Initialize output data
data_out = zeros(size(data_in, 1), L_max);

% Calculate integrated reflectance for each grating at each laser wavelength
new_R_s = cell(1, N_s);
for i = 1:N_s
    for j = 1:size(data_in, 1)
        % Integrate grating reflectance over laser spectrum
        new_R_s{i}(j) = trapz(simulation_range, R_simulated{i} .* lasers_spectra(j, :));
    end
end

% Process each wavelength
for i = 1:size(data_in, 1)
    % Process each grating
    for j = 1:N_s
        % Calculate reflection coefficient
        % For first grating, use direct reflectance
        if j == 1
            R_coef = R_simulated{j};
        else
            % For subsequent gratings, account for transmission through previous gratings
            R_coef = ones(size(R_simulated{1}));
            for k = 1:j-1
                R_coef = R_coef .* (1 - R_simulated{k});
            end
            R_coef = R_coef.^2 .* R_simulated{j};
        end
        
        % Account for crosstalk effects
        R_cx = crossOld(new_R_s{j}(i), j);
        
        % Integrate over laser spectrum
        R_coef_integrated = trapz(simulation_range, R_coef .* lasers_spectra(i, :));
        R_coef_integrated = R_coef_integrated + R_cx;
        
        % Apply fiber attenuation and grating reflection
        data_temp = R_coef_integrated * data_in(i, :) * fiber.alpha^(-2 * fbg.D(j));
        
        % Add delayed reflection to output
        data_out(i, :) = data_out(i, :) + addDelay(data_temp, D_s(j), L_max);
    end
end
end

function [data_out] = addDelay(data_in, delay, L_max)
% Add delay to signal data

data_out = zeros(size(data_in, 1), L_max);

for i = 1:size(data_in, 1)
    if L_max - delay - size(data_in, 2) < 0
        data_out(i, :) = [zeros(1, delay) data_in(i, :)];
    else
        data_out(i, :) = [zeros(1, delay) data_in(i, :) zeros(1, L_max - delay - size(data_in, 2))];
    end
end
end

function out = crossOld(R, N)
% Calculate crosstalk for a given reflectance and grating number
% This is a simplified model - actual implementation would be more complex

% Transmission coefficient
T = 1 - R;

% Simple crosstalk model based on number of gratings and reflectance
out = (N - 1) * (N - 2) / 2 * R^3 * T^(2 * N - 4);
end
