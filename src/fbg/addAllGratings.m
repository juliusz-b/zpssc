function [data_out] = addAllGratings(data_in, N_s, D_s, R_s, L_max, R_simulated, lasers_spectra, simulation_range, fbg, fiber)
% ADDALLGRATINGS Add FBG array reflections to input data, including multi-bounce ghosts
%
% [data_out] = addAllGratings(data_in, N_s, D_s, R_s, L_max, R_simulated,
%                            lasers_spectra, simulation_range, fbg, fiber)
%
% Inputs:
%   data_in          - Input optical signal data [U x N_samples]
%   N_s              - Number of gratings
%   D_s              - Grating positions (in samples)
%   R_s              - Cell array with reflectance spectrum of each grating
%   L_max            - Maximum output length
%   R_simulated      - Cell array with high-resolution reflectance spectra
%   lasers_spectra   - Laser spectral distribution [U x N_sim]
%   simulation_range - Wavelength range for simulation [1 x N_sim]
%   fbg              - FBG parameters structure
%       fbg.D            - Physical distances [m]
%       fbg.max_bounces  - (optional) Max number of reflections: 1, 3 (default 1)
%   fiber            - Fiber parameters structure
%       fiber.alpha      - Fiber attenuation (linear)
%
% Outputs:
%   data_out - Output optical signal with FBG reflections
%
% Multi-bounce model:
%   max_bounces = 1: Only primary reflections (light -> FBG_j -> detector)
%   max_bounces = 3: Primary + 3rd order ghosts (light -> FBG_j -> FBG_k -> FBG_m -> detector)
%       Ghost apparent position = D(j) - D(k) + D(m), for all j, k<j, m>k
%       Ghost amplitude = R(j)*R(k)*R(m) * T^2(through other gratings on path)

% Check input parameters
if length(R_s) < N_s
    error('Insufficient reflectance spectra provided');
end
if length(D_s) < N_s
    error('Insufficient grating positions provided');
end

% Get max_bounces parameter
if isfield(fbg, 'max_bounces')
    max_bounces = fbg.max_bounces;
else
    max_bounces = 1; % Default: primary only (backward compatible)
end

% Initialize output data
data_out = zeros(size(data_in, 1), L_max);

% Precompute: integrated reflectance for each grating at each laser wavelength
% new_R_s{j}(i) = integral of R_simulated{j} * laser_spectrum_i over wavelength
new_R_s = zeros(N_s, size(data_in, 1));
for j = 1:N_s
    for i = 1:size(data_in, 1)
        new_R_s(j, i) = trapz(simulation_range, R_simulated{j} .* lasers_spectra(i, :));
    end
end

% Precompute: cumulative transmission T^2 through gratings 1..k-1
% T2_cumul(j, :) = product of (1-R_simulated{k})^2 for k=1..j-1 (wavelength-resolved)
T2_cumul = cell(1, N_s);
T2_cumul{1} = ones(size(R_simulated{1})); % No transmission loss for first grating
for j = 2:N_s
    T2_cumul{j} = T2_cumul{j-1} .* (1 - R_simulated{j-1}).^2;
end

% ===== PRIMARY REFLECTIONS (1 bounce) =====
for i = 1:size(data_in, 1)
    for j = 1:N_s
        % Reflection coefficient: T^2(through 1..j-1) * R(j)
        R_coef = T2_cumul{j} .* R_simulated{j};

        % Integrate over laser spectrum
        R_coef_integrated = trapz(simulation_range, R_coef .* lasers_spectra(i, :));

        % Apply fiber attenuation (round-trip to grating j)
        data_temp = R_coef_integrated * data_in(i, :) * fiber.alpha^(-2 * fbg.D(j));

        % Add delayed reflection to output
        data_out(i, :) = data_out(i, :) + addDelay(data_temp, D_s(j), L_max);
    end
end

% ===== 3rd ORDER GHOSTS (3 bounces) =====
if max_bounces >= 3
    % Path: light -> FBG_j (reflect) -> FBG_k (reflect, k<j) -> FBG_m (reflect, m>k) -> detector
    % Apparent round-trip delay in samples: D_s(j) - D_s(k) + D_s(m)
    %   (because: forward D(j), back D(j)-D(k), forward D(m)-D(k), back D(m))
    %   total one-way apparent = D(j) - D(k) + D(m)
    % Amplitude: R(j) * R(k) * R(m) * transmission through other gratings on path

    neff_fiber = 1.447;
    V_fiber = 3e8 / neff_fiber;

    for i = 1:size(data_in, 1)
        for j = 1:N_s
            for k = 1:j-1  % k < j (back-reflection)
                for m = k+1:N_s  % m > k (forward re-reflection)
                    if j == m && k == j
                        continue; % Skip degenerate case
                    end

                    % Ghost apparent position in samples
                    ghost_delay = D_s(j) - D_s(k) + D_s(m);

                    % Check if ghost position fits in output
                    if ghost_delay < 0 || ghost_delay >= L_max
                        continue;
                    end

                    % Ghost amplitude (wavelength-resolved):
                    % R(j) * R(k) * R(m) * T^2 through gratings on the path
                    % Transmission: we need to account for passage through
                    % gratings between the bouncing gratings
                    R_ghost = R_simulated{j} .* R_simulated{k} .* R_simulated{m};

                    % Transmission on forward path to j: through 1..j-1
                    % Transmission on return from j to k: through k+1..j-1
                    % Transmission on forward from k to m: through k+1..m-1
                    % Transmission on return from m: through 1..m-1
                    % Total transmission factor (each passage is one-way, squared for both directions):
                    % T_fwd_to_j = prod(1-R(1..j-1))
                    % T_j_to_k   = prod(1-R(k+1..j-1))
                    % T_k_to_m   = prod(1-R(k+1..m-1))
                    % T_m_to_det = prod(1-R(1..m-1))
                    % Total = T_fwd_to_j * T_j_to_k * T_k_to_m * T_m_to_det

                    T_path = ones(size(R_simulated{1}));

                    % Forward to j: through 1..j-1
                    for g = 1:j-1
                        T_path = T_path .* (1 - R_simulated{g});
                    end

                    % j back to k: through k+1..j-1
                    for g = k+1:j-1
                        T_path = T_path .* (1 - R_simulated{g});
                    end

                    % k forward to m: through k+1..m-1
                    for g = k+1:m-1
                        T_path = T_path .* (1 - R_simulated{g});
                    end

                    % m back to detector: through 1..m-1
                    for g = 1:m-1
                        T_path = T_path .* (1 - R_simulated{g});
                    end

                    R_ghost_total = R_ghost .* T_path;

                    % Integrate over laser spectrum
                    R_ghost_integrated = trapz(simulation_range, R_ghost_total .* lasers_spectra(i, :));

                    % Fiber attenuation for the total path
                    % Total physical distance: D(j) + (D(j)-D(k)) + (D(m)-D(k)) + D(m) = 2*D(j) - 2*D(k) + 2*D(m)
                    total_dist = 2*fbg.D(j) - 2*fbg.D(k) + 2*fbg.D(m);
                    atten = fiber.alpha^(-total_dist);

                    data_temp = R_ghost_integrated * data_in(i, :) * atten;
                    data_out(i, :) = data_out(i, :) + addDelay(data_temp, ghost_delay, L_max);
                end
            end
        end
    end
end

end

function [data_out] = addDelay(data_in, delay, L_max)
% Add delay to signal data

data_out = zeros(size(data_in, 1), L_max);

for i = 1:size(data_in, 1)
    sig_len = size(data_in, 2);
    available = L_max - delay;
    if available <= 0
        % Signal doesn't fit - skip
        continue;
    end
    n_copy = min(sig_len, available);
    data_out(i, delay+1:delay+n_copy) = data_in(i, 1:n_copy);
end
end
