
%REFERENCE:A. Grabowski, J. Gustavsson, Z. S. He, and A. Larsson, “Large-Signal Equivalent Circuit for Datacom VCSELs,” Journal of Lightwave Technology,
% vol. 39, no. 10. Institute of Electrical and Electronics Engineers (IEEE), pp. 3225–3236, May 15, 2021. doi: 10.1109/jlt.2021.3064465.




% Constants
q = 1.6e-19;       % Elementary charge (C)
k = 1.38e-23;      % Boltzmann constant (J/K)
h = 6.626e-34;      % Planck constant (J·s)
c = 3e8;           % Speed of light (m/s)
lambda = 850e-9;   % Wavelength (m)
RT = 300;          % Reference temperature (K)
m0 = 9.1093837139e-31; % Electron mass (kg)

% Parameters
params.q = q;
params.k = k;
params.h = h;
params.c = c;
params.RT = RT;
params.lambda = lambda;
params.TAMB = RT; %ambient temp.

% Additional parameters
params.Aaperture = 38.5e-12;    % Aperture area (m²)
params.dA = 5 * 4e-9;           % Thickness of QW (m)
params.dB = 74e-9;              % Thickness of barrier (m)
params.dQW = 4e-9;              % QW thickness (m)
params.m0 = m0;     % Electron mass (kg)

% Thermal resistance parameters
params.rm = 170;
params.drm = -0.96;
params.drm2 = 1.38e-2;
params.drm3 = -1.11e-4;
params.Rp3 = 15;

% Bandgap energies and temperature coefficients
params.EgBRT = 1.89;
params.dEgBRTdT = -4.5e-4;
params.meB = 0.0977;          % Effective mass (relative to m0)
params.mhhB = 0.596;
params.eta_injxy = 0.8;
params.tau_sp_B = 5e-9;
params.etaB = 1;
params.eta_leak = 2;
params.A0B = 7.4e5;
params.EgCRT = 2.31;
params.dEgCRTdT = 4.94e-4;

% Active region parameters
params.EgART = 1.53;
params.dEgARTdT = -1.63e-3;
params.meA = 0.0622;
params.mhhA = 0.346;
params.tau_cap = 7e-12;
params.etaA = 1.6;
params.eta_esc = 2;

% Recombination coefficients
params.A0A = 7.4e5;
params.ART = 0.45e9;
params.dAdT = 0;
params.BRT = 3e-16;
params.dBdT = -1.19e-18;
params.dBdT2 = 5.74e-21;
params.CRT = -7.9e-42;
params.dCdT = 8.2e-44;
params.dCdT2 = -5.4e-46;

% Gain parameters
params.a0RT = 7.54e5;
params.da0dT = 8.22e2;
params.da0dT2 = 2.66;
params.a1_1 = -1.4e13;
params.a2_1 = 1.57e20;
params.a1_2 = -1.34e13;
params.a2_2 = 2.03e21;

% Carrier density parameters
params.N0RT = 2.76e24;
params.N0d1 = 9.22e21;
params.N_S = -1.5e21;
params.epsilon = 5e-24;

% Wavelength parameters
params.lambdagRT = 857e-9;
params.dlambdagdT = 0.314e-9;
params.dlambdagdnA = -3.26e-33;
params.dlambdagdnA2 = -5.44e-60;
params.lambdaRT = 847.3e-9;
params.dlambdadT = 0.064e-9;
params.v_g = 1e8;
params.Gamma = 3.28e-2;
params.beta = 4.68e-2;

% Absorption coefficients
params.alfa_absRT = 6.19e10;
params.alfa_absRTdT = 0;
params.alfa_absRTdT2 = 5e6;
params.alfa_abmRT = 1.35e10;
params.alfa_abmRTdT = -3.1e7;
params.alfa_atmRT = 1.91e11;
params.alfa_atmRTdT = -2.66e8;

% Thermal parameters
params.tau_th = 1e-7;
params.rthRT = 2.8e3;
params.drthRTdT = 2.5;