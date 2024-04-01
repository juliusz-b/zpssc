function [new_lmb, R_filtered, lmbB] = upsampleAndDenoiseGratings(R, fbg, lmb)
% UPSAMPLEANDDENOISEGRATINGS Improve FBG reflectance spectra
%
% [new_lmb, R_filtered, lmbB] = upsampleAndDenoiseGratings(R, fbg, lmb)
%
% Inputs:
%   R - FBG reflectance spectra
%   fbg - FBG parameters structure
%   lmb - Wavelength array [m]
%
% Outputs:
%   new_lmb - Upsampled wavelength array [m]
%   R_filtered - Filtered and upsampled reflectance spectra
%   lmbB - Bragg wavelengths [m]

% Set upsampling ratio
UR = 4;

% Set filter parameters
span = 20;
beta = 0.1;

% Design FIR lowpass filter
flt = designfilt('lowpassfir', 'FilterOrder', span*UR, 'CutoffFrequency', 9/33);

% Generate upsampled wavelength array
new_lmb = lmb(1):(lmb(2)-lmb(1))/UR:lmb(end)+(lmb(2)-lmb(1));

% Initialize output arrays
R_filtered = [];
lmbB = [];

% Process each grating reflectance spectrum
for i = 1:size(R, 1)
    % Get current reflectance spectrum
    R_temp = R(i, :);
    
    % Upsample the spectrum
    R_upsampled = upsample(R_temp, UR);
    
    % Apply filtering to smooth the upsampled spectrum
    R_filtered(i, :) = abs(filterMySignal(R_upsampled, flt.Coefficients) * UR);
    
    % Calculate Bragg wavelength as weighted average
    lmbB(i) = sum(R_filtered(i, :) .* new_lmb) / sum(R_filtered(i, :));
end
end
