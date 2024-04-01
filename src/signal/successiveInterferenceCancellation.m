function [data_out, detected_signals] = successiveInterferenceCancellation(data_in, reference_data, N_iterations)
% SUCCESSIVEINTERFERENCECANCELLATION Implement SIC algorithm for FBG detection
%
% [data_out, detected_signals] = successiveInterferenceCancellation(data_in, reference_data, N_iterations)
%
% Implements the Successive Interference Cancellation (SIC) algorithm
% to improve detection in the presence of multiple interfering signals.
%
% Inputs:
%   data_in - Input signal data
%   reference_data - Reference signals (spreading codes)
%   N_iterations - Number of SIC iterations
%
% Outputs:
%   data_out - Processed signal after interference cancellation
%   detected_signals - Information about detected signals

% Initialize output
data_out = data_in;
detected_signals = struct('position', {}, 'amplitude', {}, 'code_index', {});

% Perform SIC iterations
for iter = 1:N_iterations
    % Calculate cross-correlation for all reference signals
    xcorr_data = zeros(size(reference_data, 1), length(data_out));
    for i = 1:size(reference_data, 1)
        xcorr_temp = xcorr(data_out, reference_data(i, :));
        xcorr_data(i, :) = xcorr_temp(ceil(length(xcorr_temp)/2):end);
    end
    
    % Find strongest signal
    [max_vals, max_locs] = max(xcorr_data, [], 2);
    [strongest_val, strongest_code] = max(max_vals);
    strongest_loc = max_locs(strongest_code);
    
    % If the strongest signal is below threshold, stop iterations
    threshold = 0.1 * max(data_in);
    if strongest_val < threshold
        break;
    end
    
    % Store detected signal information
    detected_signals(iter).position = strongest_loc;
    detected_signals(iter).amplitude = strongest_val;
    detected_signals(iter).code_index = strongest_code;
    
    % Reconstruct the detected signal
    detected_signal = zeros(size(data_out));
    if strongest_loc + length(reference_data(strongest_code, :)) - 1 <= length(detected_signal)
        detected_signal(strongest_loc:strongest_loc + length(reference_data(strongest_code, :)) - 1) = ...
            strongest_val * reference_data(strongest_code, :);
    else
        overlap = strongest_loc + length(reference_data(strongest_code, :)) - 1 - length(detected_signal);
        detected_signal(strongest_loc:end) = ...
            strongest_val * reference_data(strongest_code, 1:end-overlap);
    end
    
    % Subtract the detected signal from the input
    data_out = data_out - detected_signal;
    
    % Calculate SINR improvement
    SINR_before = 10 * log10(strongest_val^2 / var(data_in - detected_signal));
    SINR_after = 10 * log10(strongest_val^2 / var(data_out));
    
    fprintf('SIC Iteration %d: Removed signal at position %d with amplitude %.4f (SINR improvement: %.2f dB)\n', ...
        iter, strongest_loc, strongest_val, SINR_after - SINR_before);
end

end
