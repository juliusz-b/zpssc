function [data_out] = noiseAndMeanDataInChannel(data_in, K, laser, noise, pd)
% NOISEANDMEANDATAINCHANNEL Add noise and average multiple measurements
%
% [data_out] = noiseAndMeanDataInChannel(data_in, SNR, K, laser, noise, pd)
%
% Inputs:
%   data_in - Input signal data
%   SNR - Signal-to-noise ratio [dB]
%   K - Number of averages
%   laser - Laser parameters structure
%   noise - Noise parameters structure
%   pd - Photodetector parameters structure
%
% Outputs:
%   data_out - Averaged signal with added noise

% Initialize output with zeros
data_out = zeros(size(data_in));

% Add noise and accumulate for K iterations
for i = 1:K
    data_out = data_out + addNoise(data_in, laser, noise, pd);
end

% Average the accumulated data
data_out = data_out / K;
end
