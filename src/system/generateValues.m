function [D_s, L_max, time_domain, length_domain] = generateValues(N, Fsample, Nb, D_s, U, N_s, should_add_extra_zeros)
% GENERATEVALUES Generate parameter values for simulation
%
% [D_s, L_max, time_domain, length_domain] = generateValues(N, Fsample, Nb, D_s, U, N_s, should_add_extra_zeros)
%
% Inputs:
%   N - Code length
%   Fsample - Sampling frequency [Hz]
%   Nb - Samples per bit
%   D_s - FBG distances [m]
%   U - Number of wavelength channels
%   N_s - Number of gratings
%   should_add_extra_zeros - (optional) Flag to add extra zeros padding
%
% Outputs:
%   D_s - FBG positions in samples
%   L_max - Maximum signal length
%   time_domain - Time domain array [s]
%   length_domain - Length domain array [m]

% Set default for extra zeros
if nargin < 7
   should_add_extra_zeros = false; 
end

% Double sample count if extra zeros requested
if should_add_extra_zeros
    Nb = Nb * 2;
end

% Check if all FBG distances are provided
if length(D_s) < N_s
    error('FBG distances not fully specified');
end

% Calculate fiber parameters
neff = 1.447;  % Effective refractive index
V = 3e8 / neff; % Light velocity in fiber

% Calculate measurement frequency
f_measure = Fsample / N / Nb / 32;

% Display key parameters
disp(['Minimum code length required: ' sprintf('%i', ceil(4 * (max([N_s .* D_s])) * (Fsample / Nb) / V)) ' samples']);
disp(['No signal overlap when gratings separated by at least ' num2str(ceil(Nb * N / Fsample * V / 2)) ' m']);

% Convert physical distances to sample counts
D_s = ceil(2 * D_s / V * Fsample);

% Determine maximum signal length needed
L_max = Nb * N + ceil(max(D_s(1:N_s)));

% Ensure even length for better FFT performance
if mod(L_max, 2) == 1
    L_max = L_max + 1;
end

% Generate time and length domains
time_domain = (0:L_max-1) / Fsample;
length_domain = (0:L_max-1) ./ (2 / V * Fsample);

% Display system resolution information
disp(['Spatial resolution: ' sprintf('%.2f m', V / Fsample / 2)]);
disp(['Sampling frequency: ' sprintf('%.2f MHz', Fsample * 1e-6)]);
disp(['Measurement frequency: ' sprintf('%.2f kHz', f_measure * 1e-3)]);
end
