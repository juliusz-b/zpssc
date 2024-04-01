function [R, T] = generateBraggGrating(lambdas, period, neff, deltaneff, grating_length, visibility, sections)
% GENERATEBRAGGGRATING Simulates the spectral response of a fiber Bragg grating
%
% [R, T] = generateBraggGrating(lambdas, period, neff, deltaneff, grating_length, visibility, sections)
%
% Inputs:
%   lambdas - Wavelength array for simulation [m]
%   period - Grating period [m]
%   neff - Effective refractive index
%   deltaneff - Refractive index modulation depth
%   grating_length - Length of the grating [m]
%   visibility - Fringe visibility (typically 1)
%   sections - Number of sections to divide the grating for simulation
%
% Outputs:
%   R - Reflectance spectrum
%   T - Transfer matrix for each wavelength

% Calculate Bragg wavelength
lambda_bragg = 2 * period * neff;

% Calculate length of each section
L = grating_length / sections;

% Initialize transfer matrix
T = repmat(eye(2), 1, length(lambdas));

% Generate apodization profile (Gaussian)
apod = gaussmf(1:sections, [sections/3 sections/2]);

% Calculate transfer matrix for each section
for i = 1:sections
    % Calculate detuning parameter
    sigma = 2 * pi ./ lambdas * deltaneff * apod(i) + ...
            2 * pi * neff * (1 ./ lambdas - 1 / lambda_bragg);
    
    % Calculate coupling coefficient
    kappa = pi * visibility * deltaneff * apod(i) ./ lambdas;
    
    % Calculate gamma parameter
    gamma = sqrt(kappa.^2 - sigma.^2);
    
    % Calculate transfer matrix components
    T11 = cosh(gamma * L) - 1i * sigma ./ gamma .* sinh(gamma * L);
    T12 = -1i * kappa ./ gamma .* sinh(gamma * L);
    T21 = 1i * kappa ./ gamma .* sinh(gamma * L);
    T22 = cosh(gamma * L) + 1i * sigma ./ gamma .* sinh(gamma * L);
    
    % Combine matrices
    T_temp = matrixCombine(T11, T12, T21, T22);
    T = matrixMult(T_temp, T, 2);
end

% Calculate reflectance from transfer matrix
R = 1 - abs(1 ./ T(1, 1:2:size(T, 2))).^2;
end

function [C] = matrixCombine(A, B, C, D)
% Combine matrix components into block matrices
if ~isequal(size(A), size(B), size(C), size(D))
    error('Wrong dimensions!');
end
if size(A, 1) > size(A, 2)
    A = A';
    B = B';
    C = C';
    D = D';
end

E = zeros(2, size(A, 2) * 2);

for i = 1:size(A, 2)
    E(:, (i-1)*2+1:(i-1)*2+2) = [A(i), B(i); C(i), D(i)];
end
C = E;
end

function [C] = matrixMult(A, B, k)
% Multiply block matrices
if nargin < 3
   k = 2; 
end

if length(size(A)) ~= 2 || length(size(B)) ~= 2
   error('Wrong dimensions!'); 
end

if sum(size(A) == size(B)) ~= length(size(A))
   error('Wrong dimensions!'); 
end

if size(A, 1) > size(A, 2)
    A = A';
end
if size(B, 1) > size(B, 2)
    B = B';
end
    
N = size(A, 2) / k;
if N ~= floor(N)
   error('Wrong dimensions!'); 
end

C = zeros(size(A));

for i = 1:N
    C(:, (i-1)*k+1:(i-1)*k+k) = A(:, (i-1)*k+1:(i-1)*k+k) * B(:, (i-1)*k+1:(i-1)*k+k);
end
end
