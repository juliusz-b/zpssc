function [data_out] = addNoise(data_in, snr, laser, noise, pd)
% ADDNOISE Add realistic noise to optical signal
%
% [data_out] = addNoise(data_in, snr, laser, noise, pd)
%
% Inputs:
%   data_in - Input signal data
%   snr - Signal-to-noise ratio [dB]
%   laser - Laser parameters structure
%   noise - Noise parameters structure:
%       noise.type - 'awgn-snr', 'snr-relative', 'nep', 'true'
%       noise.SNR - SNR for 'awgn-snr' and 'snr-relative'
%       noise.NEP - Noise equivalent power [W/√Hz] for 'nep'
%       noise.Fn - Noise figure for 'true'
%   pd - Photodetector parameters structure:
%       pd.A - Responsivity [A/W]
%       pd.BW - Bandwidth [Hz]
%       pd.RL - Load resistance [Ω]
%       pd.Idark - Dark current [A]
%
% Outputs:
%   data_out - Signal with added noise

% Initialize output
data_out = data_in;

% Check for complex