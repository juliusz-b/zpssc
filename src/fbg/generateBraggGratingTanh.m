function [R] = generateBraggGratingTanh(lambdas, period, neff, deltaneff, grating_length, visibility)
% GENERATEBRAGGGRATINGTANH Analytical tanh model for uniform FBG reflectance
%
% [R] = generateBraggGratingTanh(lambdas, period, neff, deltaneff, grating_length, visibility)
%
% Inputs:
%   lambdas        - Wavelength array [m]
%   period         - Grating period [m]
%   neff           - Effective refractive index
%   deltaneff      - Refractive index modulation depth
%   grating_length - Length of the grating [m]
%   visibility     - Fringe visibility (typically 1)
%
% Outputs:
%   R - Reflectance spectrum
%
% Uses coupled-mode theory analytical solution for uniform FBG:
%   R = tanh^2(kappa * L) at Bragg wavelength (peak)
%   General: R(lambda) = kappa^2 * sinh^2(gamma*L) / (gamma^2 * cosh^2(gamma*L) + sigma^2 * sinh^2(gamma*L))
%   where gamma = sqrt(kappa^2 - sigma^2), sigma = detuning

% Bragg wavelength
lambda_bragg = 2 * neff * period;

% Coupling coefficient (wavelength-dependent)
kappa = pi * visibility * deltaneff ./ lambdas;

% Detuning parameter (no apodization - uniform grating)
sigma = 2 * pi * neff * (1 ./ lambdas - 1 / lambda_bragg);

% Propagation constant
gamma2 = kappa.^2 - sigma.^2;

% Initialize reflectance
R = zeros(size(lambdas));

% Case 1: kappa^2 > sigma^2 (near Bragg wavelength) - real gamma
idx_real = gamma2 > 0;
gamma_r = sqrt(gamma2(idx_real));
L = grating_length;
sinh_gL = sinh(gamma_r * L);
cosh_gL = cosh(gamma_r * L);
R(idx_real) = kappa(idx_real).^2 .* sinh_gL.^2 ./ ...
    (gamma_r.^2 .* cosh_gL.^2 + sigma(idx_real).^2 .* sinh_gL.^2);

% Case 2: kappa^2 < sigma^2 (far from Bragg) - imaginary gamma -> sin/cos
idx_imag = gamma2 < 0;
gamma_i = sqrt(-gamma2(idx_imag));
sin_gL = sin(gamma_i * L);
cos_gL = cos(gamma_i * L);
R(idx_imag) = kappa(idx_imag).^2 .* sin_gL.^2 ./ ...
    (gamma_i.^2 .* cos_gL.^2 + sigma(idx_imag).^2 .* sin_gL.^2);

% Case 3: exactly at boundary (gamma = 0) - L'Hôpital
idx_zero = gamma2 == 0;
R(idx_zero) = kappa(idx_zero).^2 * L^2 ./ (1 + sigma(idx_zero).^2 * L^2);

% Clamp to [0, 1]
R = max(0, min(1, R));

end
