% Large-Signal VCSEL Model with Scaled Variables

% Clear workspace and figures
clear; clc; close all;

% Constants
q = 1.602176634e-19;     % Elementary charge (C)
k = 1.380649e-23;        % Boltzmann constant (J/K)
h = 6.62607015e-34;      % Planck constant (J·s)
c = 2.99792458e8;        % Speed of light (m/s)
lambda = 850e-9;         % Wavelength (m)
RT = 300;                % Reference temperature (K)
Ibias = 1e-3;            % Bias current (A)

% Scaling Factors
N0 = 1e24;               % Reference carrier concentration (m^-3)
S0 = 1e12;               % Reference photon number
T0 = RT;                 % Reference temperature (K)
V0 = 1e-18;              % Reference volume (m^3)
E0 = q;                  % Reference energy (J)
t0 = 1e-9;               % Reference time (s)

% Parameters (scaled where appropriate)
params.q = q;
params.k = k;
params.h = h;
params.c = c;
params.RT = RT / T0;                   % Scaled reference temperature
params.lambda = lambda;
params.N0 = N0;
params.S0 = S0;
params.T0 = T0;

% Additional parameters (scaled)
params.Aaperture = 38.5e-12;           % Aperture area (m²)
params.dA = 5 * 4e-9;                  % Thickness of QW (m)
params.dB = 74e-9;                     % Thickness of barrier (m)
params.dQW = 4e-9;                     % QW thickness (m)
params.VA = params.dA * params.Aaperture / V0;  % Scaled volume of active region
params.VB = params.dB * params.Aaperture / V0;  % Scaled volume of barrier region
params.m0 = 9.10938356e-31;            % Electron rest mass (kg)

% Thermal resistance parameters
params.rm = 170;
params.drm = -0.96;
params.drm2 = 1.38e-2;
params.drm3 = -1.11e-4;

% Bandgap energies and temperature coefficients (in eV)
params.EgBRT = 1.89;
params.dEgBRTdT = -4.5e-4;
params.meB = 0.0977 * params.m0;       % Effective electron mass in barriers (kg)
params.mhhB = 0.596 * params.m0;       % Effective hole mass in barriers (kg)
params.eta_injxy = 0.8;
params.tau_sp_B = 5e-9 / t0;           % Scaled spontaneous recombination lifetime
params.etaB = 1;
params.eta_leak = 2;
params.A0B = 7.4e5;
params.EgCRT = 2.31;
params.dEgCRTdT = 4.94e-4;

% Active region parameters
params.EgART = 1.53;
params.dEgARTdT = -1.63e-3;
params.meA = 0.0622 * params.m0;       % Effective electron mass in QWs (kg)
params.mhhA = 0.346 * params.m0;       % Effective hole mass in QWs (kg)
params.tau_cap = 7e-12 / t0;           % Scaled capture time
params.etaA = 1.6;
params.eta_esc = 2;

% Recombination coefficients (scaled)
params.A0A = 7.4e5;
params.ART = 0.45e9 * t0;             % Scaled A coefficient
params.dAdT = 0;
params.BRT = 3e-16 * t0 * N0;         % Scaled B coefficient
params.dBdT = -1.19e-18 * t0 * N0;
params.dBdT2 = 5.74e-21 * t0 * N0;
params.CRT = -7.9e-42 * t0 * N0^2;    % Scaled C coefficient
params.dCdT = 8.2e-44 * t0 * N0^2;
params.dCdT2 = -5.4e-46 * t0 * N0^2;

% Gain parameters (scaled)
params.a0RT = 7.54e5 / N0;            % Scaled gain coefficient
params.da0dT = 8.22e2 / N0;
params.da0dT2 = 2.66 / N0;
params.a1_1 = -1.4e13 / N0;
params.a2_1 = 1.57e20 / N0^2;
params.a1_2 = -1.34e13 / N0;
params.a2_2 = 2.03e21 / N0^2;

% Carrier density parameters
params.N0RT = 2.76e24 / N0;           % Scaled transparency carrier density
params.N0d1 = 9.22e21 / N0;
params.N_S = -1.5e21 / N0;
params.epsilon = 5e-24 * N0 / S0;

% Wavelength parameters
params.lambdagRT = 857e-9;
params.dlambdagdT = 0.314e-9;
params.dlambdagdnA = -3.26e-33 * N0;
params.dlambdagdnA2 = -5.44e-60 * N0^2;
params.lambdaRT = 847.3e-9;
params.dlambdadT = 0.064e-9;
params.v_g = 1e8 / c;                 % Scaled group velocity
params.Gamma = 3.28e-2;
params.beta = 4.68e-2;

% Absorption coefficients (scaled)
params.alfa_absRT = 6.19e10 * t0;
params.alfa_absRTdT = 0;
params.alfa_absRTdT2 = 5e6 * t0;
params.alfa_abmRT = 1.35e10 * t0;
params.alfa_abmRTdT = -3.1e7 * t0;
params.alfa_atmRT = 1.91e11 * t0;
params.alfa_atmRTdT = -2.66e8 * t0;

% Thermal parameters
params.tau_th = 1e-7 / t0;            % Scaled thermal time constant
params.rthRT = 2.8e3;
params.drthRTdT = 2.5;

% Initial Conditions (scaled)
nB0 = (Ibias / q / params.eta_injxy) / N0;  % Scaled barrier carrier concentration
nA0 = params.N0RT;                          % Scaled QW carrier concentration
initial_conditions = [nB0, nA0, 1e-3 / S0, 0];  % [nB_hat, nA_hat, S_hat, dT_hat]

% Driving current signal (scaled)
time_vector = linspace(0, 1e-9, 1000) / t0;     % Scaled time vector
current_signal = Ibias + (Ibias / 2) * sin(2 * pi * 1e9 * time_vector * t0);  % Signal
current_signal = current_signal / (q * N0 * V0 / t0);  % Scaled current

% Interpolated current function
current_func = @(t) interp1(time_vector, current_signal, t, 'linear', 'extrap');

% Solve the ODEs
tspan = [time_vector(1), time_vector(end)];  % Simulation time span
options = odeset('RelTol',1e-6,'AbsTol',1e-9,'MaxStep',1e-2);
[t, y] = ode45(@(t, y) vcsel_odes(t, y, params, current_func), tspan, initial_conditions, options);

% Extract results (unscale)
nB = y(:, 1) * N0;  % Carrier concentration in barriers
nA = y(:, 2) * N0;  % Carrier concentration in QWs
S = y(:, 3) * S0;   % Photon number
T = y(:, 4) * T0 + RT;  % Temperature in K

% Calculate optical power (placeholder calculation)
P_opt = S;  % Adjust according to actual optical power calculation

% Plot results
figure;
subplot(3, 1, 1);
plot(t * t0, P_opt);
xlabel('Time (s)');
ylabel('Optical Power (arbitrary units)');
title('Optical Power');

subplot(3, 1, 2);
plot(t * t0, nA);
xlabel('Time (s)');
ylabel('Carrier Concentration in QWs (m^{-3})');
title('Carrier Concentration in QWs');

subplot(3, 1, 3);
plot(t * t0, T);
xlabel('Time (s)');
ylabel('Temperature (K)');
title('Internal Temperature');

% Function definitions
function dydt = vcsel_odes(t, y, params, current_func)
    % Extract scaled variables
    nB_hat = y(1);  % Scaled barrier carrier concentration
    nA_hat = y(2);  % Scaled QW carrier concentration
    S_hat = y(3);   % Scaled photon number
    dT_hat = y(4);  % Scaled temperature difference
    
    % Unscale variables
    N0 = params.N0;
    S0 = params.S0;
    T0 = params.T0;
    V0 = 1e-18;    % Reference volume (m^3)
    t0 = 1e-9;     % Reference time (s)
    
    nB = nB_hat * N0;
    nA = nA_hat * N0;
    S = S_hat * S0;
    T = (dT_hat * T0) + params.RT * T0;
    
    % Extract parameters
    c = params.c;
    q = params.q;
    k = params.k;
    h = params.h;
    RT = params.RT * T0;
    tau_sp_B = params.tau_sp_B * t0;
    tau_cap = params.tau_cap * t0;
    VA = params.VA * V0;
    VB = params.VB * V0;
    N_S = params.N_S * N0;
    epsilon = params.epsilon * S0 / N0;
    
    % Effective density of states
    NCT = 2 * (2 * pi * params.meB * k * T / h^2)^(3/2);
    NVT = 2 * (2 * pi * params.mhhB * k * T / h^2)^(3/2);
    
    % Bandgap energies
    Eg_A = calcEgA(T, RT, params.EgART, params.dEgARTdT);
    Eg_B = calcEgB(T, RT, params.EgBRT, params.dEgBRTdT);
    Eg_C = calcEgC(T, RT, params.EgCRT, params.dEgCRTdT);
    
    % Recombination coefficients
    A = calcA(T, RT, params.ART, params.dAdT);
    B = calcB(T, RT, params.BRT, params.dBdT, params.dBdT2);
    C = calcC(T, RT, params.CRT, params.dCdT, params.dCdT2);
    N0_param = calcN0(T, RT, params.N0RT * N0, params.N0d1 * N0);
    a0 = calca0(T, RT, params.a0RT * N0, params.da0dT * N0, params.da0dT2 * N0);
    
    lambda = calclambda(T, RT, params.lambdaRT, params.dlambdadT);
    lambdag = calclambdag(T, RT, nA, VA, params.lambdagRT, params.dlambdagdT, params.dlambdagdnA, params.dlambdagdnA2);
    
    % Gain parameters
    if lambda <= lambdag
        a1 = params.a1_1 * N0;
        a2 = params.a2_1 * N0^2;
    else
        a1 = params.a1_2 * N0;
        a2 = params.a2_2 * N0^2;
    end
    
    % Quasi-Fermi levels (in J)
    epsilon_small = 1e-20;  % Avoid log(0)
    EfA_numerator = h^2 * params.dQW * nA / VA;
    EfA_denominator_e = 4 * pi * k * T * params.meA;
    EfA_denominator_h = 4 * pi * k * T * params.mhhA;
    EfA_term_e = log(max(exp(EfA_numerator / EfA_denominator_e) - 1, epsilon_small));
    EfA_term_h = log(max(exp(EfA_numerator / EfA_denominator_h) - 1, epsilon_small));
    EfA = Eg_A * q + params.etaA * k * T * (EfA_term_e + EfA_term_h);
    
    % Barrier quasi-Fermi levels
    var1C = nB / (VB * NCT);
    var1V = nB / (VB * NVT);
    var2C = 3 * sqrt(pi) * nB / (VB * (4 * NCT)^(2/3));
    var2V = 3 * sqrt(pi) * nB / (VB * (4 * NVT)^(2/3));
    
    el1C = log(max(var1C, epsilon_small)) / (1 - var1C^2);
    el2C = var2C / (1 + (0.24 + 1.08 * var2C)^(-2));
    el1V = log(max(var1V, epsilon_small)) / (1 - var1V^2);
    el2V = var2V / (1 + (0.24 + 1.08 * var2V)^(-2));
    
    EfB = Eg_B * q + params.etaB * k * T * (el1C + el2C + el1V + el2V);
    
    % Absorption coefficients
    alpha_tm = calcalfaatm(T, RT, params.alfa_atmRT, params.alfa_atmRTdT);
    alpha_abs = calcalfaabs(T, RT, params.alfa_absRT, params.alfa_absRTdT, params.alfa_absRTdT2);
    alpha_bm = calcalfaabm(T, RT, params.alfa_abmRT, params.alfa_abmRTdT);
    
    eta_inj = params.eta_injxy;
    rth = calcrth(T, RT, params.rthRT, params.drthRTdT);
    tau_th = params.tau_th * t0;
    
    % Drive current and internal current
    ISCH = current_func(t);  % Interpolated current (scaled)
    i_inj = eta_inj * ISCH;
    
    % Barrier dynamics
    i_cap = nB / tau_cap;
    i_sp_B = nB / tau_sp_B;
    
    % QW dynamics
    i_esc = params.A0A * VA * T^2 / (q * params.dA) * exp(-Eg_B / (params.eta_esc * k * T / q)) * (exp(EfA / (params.eta_esc * k * T)) - 1);
    i_leak = params.A0B * VB * T^2 / (q * params.dB) * exp(-Eg_C / (params.eta_leak * k * T / q)) * (exp(EfB / (params.eta_leak * k * T)) - 1);
    if isnan(i_esc) || isinf(i_esc), i_esc = 0; end
    if isnan(i_leak) || isinf(i_leak), i_leak = 0; end
    
    i_sp_A = A * nA + B * nA^2 / VA + C * nA^3 / VA^2;
    
    % Gain
    g_num = (nA / VA + N_S) / (N0_param + N_S);
    if g_num <= 0, g_num = epsilon_small; end  % Avoid log of non-positive
    g = (a0 * log(g_num) - a1 * (lambda - lambdag) - a2 * (lambda - lambdag)^2) * ...
        (1 / (1 + epsilon * S * params.Gamma / VA));
    i_st = params.Gamma * params.v_g * g * S;
    
    % Photon dynamics
    i_abs = alpha_abs * S;
    i_tm = alpha_tm * S;
    isp = params.beta * B * nA^2 / VA;
    i_bm = alpha_bm * S;
    
    % Energy calculations (simplified)
    USCH = (nA * EfA + nB * EfB) / (nA + nB);
    
    % Temperature dynamics
    g_gen = ISCH * USCH - i_tm * h * c / lambda;
    g_diss = (T - RT) / rth;
    dTdt = (g_gen - g_diss) / (tau_th * rth);
    
    % Rate equations
    dnBdt = i_inj - i_sp_B - i_cap + i_esc - i_leak;
    dnAdt = i_cap - i_sp_A - i_esc - i_st;
    dSdt = isp + i_st - i_abs - i_bm - i_tm;
    
    % Scale derivatives
    dnBdt_hat = dnBdt * t0 / N0;
    dnAdt_hat = dnAdt * t0 / N0;
    dSdt_hat = dSdt * t0 / S0;
    ddTdt_hat = dTdt * t0 / T0;
    
    % Combine derivatives
    dydt = [dnBdt_hat; dnAdt_hat; dSdt_hat; ddTdt_hat]
end

% Auxiliary functions
function rm = calcrm(T, RT, rmRT, drm, drm2, drm3)
    rm = rmRT + (T - RT) * drm + 0.5 * (T - RT)^2 * drm2 + (1 / 6) * (T - RT)^3 * drm3;
end

function EgA = calcEgA(T, RT, EgART, dEgARTdT)
    EgA = EgART + (T - RT) * dEgARTdT;
end

function EgB = calcEgB(T, RT, EgBRT, dEgBRTdT)
    EgB = EgBRT + (T - RT) * dEgBRTdT;
end

function EgC = calcEgC(T, RT, EgCRT, dEgCRTdT)
    EgC = EgCRT + (T - RT) * dEgCRTdT;
end

function A = calcA(T, RT, ART, dAdT)
    A = ART + (T - RT) * dAdT;
end

function B = calcB(T, RT, BRT, dBdT, dBdT2)
    B = BRT + (T - RT) * dBdT + 0.5 * (T - RT)^2 * dBdT2;
end

function C = calcC(T, RT, CRT, dCdT, dCdT2)
    C = CRT + (T - RT) * dCdT + 0.5 * (T - RT)^2 * dCdT2;
end

function a0 = calca0(T, RT, a0RT, da0dT, da0dT2)
    a0 = a0RT + (T - RT) * da0dT + 0.5 * (T - RT)^2 * da0dT2;
end

function N0 = calcN0(T, RT, N0RT, N0d1)
    N0 = N0RT + (T - RT) * N0d1;
end

function lambda = calclambda(T, RT, lambdaRT, dlambdadT)
    lambda = lambdaRT + (T - RT) * dlambdadT;
end

function lambdag = calclambdag(T, RT, nA, VA, lambdagRT, dlambdagdT, dlambdagdnA, dlambdagdnA2)
    nA_VA = nA / VA;
    lambdag = lambdagRT + (T - RT) * dlambdagdT + nA_VA * dlambdagdnA + 0.5 * nA_VA^2 * dlambdagdnA2;
end

function alfa_abs = calcalfaabs(T, RT, alfa_absRT, alfa_absRTdT, alfa_absRTdT2)
    alfa_abs = alfa_absRT + (T - RT) * alfa_absRTdT + 0.5 * (T - RT)^2 * alfa_absRTdT2;
end

function alfa_abm = calcalfaabm(T, RT, alfa_abmRT, alfa_abmRTdT)
    alfa_abm = alfa_abmRT + (T - RT) * alfa_abmRTdT;
end

function alfa_atm = calcalfaatm(T, RT, alfa_atmRT, alfa_atmRTdT)
    alfa_atm = alfa_atmRT + (T - RT) * alfa_atmRTdT;
end

function rth = calcrth(T, RT, rthRT, drthRTdT)
    rth = rthRT + (T - RT) * drthRTdT;
end
