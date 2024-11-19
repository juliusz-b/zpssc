function [syg] = VCSELtransmission1(x_in,SR,Ibias,Imod)
% Large-Signal VCSEL Model with Signal Vector Input

%REFERENCE:A. Grabowski, J. Gustavsson, Z. S. He, and A. Larsson, “Large-Signal Equivalent Circuit for Datacom VCSELs,” Journal of Lightwave Technology,
% vol. 39, no. 10. Institute of Electrical and Electronics Engineers (IEEE), pp. 3225–3236, May 15, 2021. doi: 10.1109/jlt.2021.3064465.

%Ibias = 10e-3;      % Bias current (A)
%Imod = 5e-3;


%SR = 400e9;
%SRCAP = SR;
%BR = SRCAP/UF;
NLED = length(x_in);
time_vector = 0:1/SR:(1/SR)*(NLED-1);

%time_vector = linspace(0, 100e-9, 10000);  % Time vector in seconds
%current_signal = Ibias + (Ibias / 2) * sin(2 * pi * 1e9 * time_vector);  % Example signal
current_signal = Ibias + (Imod) * x_in;  % Example signal
%current_signal = Ibias + 0*(Imod) * ones(size(x_in));  % Example signal
%current_signal = [ones(length(x_in),1)*Ibias; zeros(length(x_in)/2,1); ones(length(x_in)/2,1)*Ibias];

VCSELparams1;

% Initial Conditions
nB0 = Ibias / q * params.eta_injxy / (1/params.tau_sp_B+1/params.tau_cap);
%nA0 = params.N0RT;  % Approximate with transparency carrier density
nA0 = (-params.ART +sqrt(params.ART^2-4*params.BRT/(params.dA * params.Aaperture)*(-nB0/params.tau_cap) ) )/(2*params.BRT/(params.dA * params.Aaperture));
S0 = params.beta*params.BRT*nA0^2/(params.dA * params.Aaperture)/(params.alfa_absRT+params.alfa_abmRT+params.alfa_atmRT-params.Gamma*params.v_g*1);
%initial_conditions = [nB0, nA0, 1e-9, RT];  % [nB, nA, S, T]
initial_conditions = [nB0 nA0 S0 RT];  % [nB, nA, S, T]
initial_conditions = [0 0 0 RT];  % [nB, nA, S, T]

% Interpolated current function
%current_func = @(t) interp1q(time_vector, current_signal, t);
current_func = @(t) interp1(time_vector, current_signal, t, 'linear', 'extrap');


tspan = [time_vector(1), time_vector(end)];
tic
[t, y] = ode15s(@(t, y) vcsel_odes(t, y, params, current_func), tspan, initial_conditions);
toc

nB = y(:, 1);  % Carrier concentration in barriers
nA = y(:, 2);  % Carrier concentration in QWs
S = y(:, 3);   % Photon number
T = y(:, 4);   % Temperature


% opt out
itm = S.*calcalfaatm(T, RT, params.alfa_atmRT, params.alfa_atmRTdT);
P_opt = itm*h*c./calclambda(T, RT, params.lambdaRT, params.dlambdadT);
P_opt = interp1(t,P_opt,time_vector);
nA = interp1(t,nA,time_vector);
nB = interp1(t,nB,time_vector);
T = interp1(t,T,time_vector);
t = time_vector;


U = calcUdiode(current_signal,T,RT,nA,nB,params);


syg = P_opt';


% Plot results
figure;
subplot(3, 2, 1);
plot(t, P_opt*1e3);
xlabel('Time (s)');
ylabel('Optical Power (mW)');
title('Optical Power');

subplot(3, 2, 2);
plot(t, nA);
xlabel('Time (s)');
ylabel('Carrier Concentration in QWs (m^{-3})');
title('Carrier Concentration in QWs');

subplot(3, 2, 3);
plot(t, T-275);
xlabel('Time (s)');
ylabel('Temperature (^oC)');
title('Internal Temperature');


subplot(3, 2, 4);
plot(t, U);
xlabel('Time (s)');
ylabel('Voltage drop (V)');
title('Voltage drop');

subplot(3, 2, 5);
plot(t, calclambda(T, RT, params.lambdaRT, params.dlambdadT)*1e9);
xlabel('Time (s)');
ylabel('Output wavelength (nm)');
title('Output wavelength');










end


% Function definitions
function dydt = vcsel_odes(t, y, params, current_func)
    % Extract variables
    nB = y(1);  % Carrier concentration in barriers
    nA = y(2);  % Carrier concentration in QWs
    S = y(3);   % Photon number
    T = y(4);   % Temperature
    
    if ~isreal(y)
        aaa = 0;
    end

    % Extract parameters
    c = params.c;
    q = params.q;
    k = params.k;
    h = params.h;
    RT = params.RT;
    tau_sp_B = params.tau_sp_B;
    tau_cap = params.tau_cap;
    VA = params.dA * params.Aaperture;  % Volume of active region (QWs)
    VB = params.dB * params.Aaperture;  % Volume of barrier region
    N_S = params.N_S;
    epsilon = params.epsilon;


    % Bandgap energies
    Eg_A = calcEgA(T, RT, params.EgART, params.dEgARTdT);
    Eg_B = calcEgB(T, RT, params.EgBRT, params.dEgBRTdT);
    Eg_C = calcEgC(T, RT, params.EgCRT, params.dEgCRTdT);

    % Recombination coefficients
    A = calcA(T, RT, params.ART, params.dAdT);
    B = calcB(T, RT, params.BRT, params.dBdT, params.dBdT2);
    C = calcC(T, RT, params.CRT, params.dCdT, params.dCdT2);
    N0 = calcN0(T, RT, params.N0RT, params.N0d1);
    a0 = calca0(T, RT, params.a0RT, params.da0dT, params.da0dT2);

    lambda = calclambda(T, RT, params.lambdaRT, params.dlambdadT);
    lambdag = calclambdag(T, RT, nA, VA, params.lambdagRT, params.dlambdagdT, params.dlambdagdnA, params.dlambdagdnA2);

    % Gain parameters
    if lambda <= lambdag
        a1 = params.a1_1;
        a2 = params.a2_1;
    else
        a1 = params.a1_2;
        a2 = params.a2_2;
    end

    % Quasi-Fermi levels
    %EfA = Eg_A + params.etaA *  ( ...
    %    k * T *log(exp(h^2 * params.dQW * nA / VA / (4 * pi * params.k * T * params.m0 * params.meA)) - 1) + ...
    %    k * T *log(exp(h^2 * params.dQW * nA / VA / (4 * pi * params.k * T * params.m0 * params.mhhA)) - 1) ...
    %);
    EfA = calcEfA(T,RT,nA,params);
    % 
    % % Define x_e and x_h
    % x_e = h^2 * params.dQW * nA / (VA * (4 * pi * k * T * params.m0 * params.meA));
    % x_h = h^2 * params.dQW * nA / (VA * (4 * pi * k * T * params.m0 * params.mhhA));
    % 
    % % Avoid numerical issues
    % epsilon_small = 1e-20;
    % x_e = max(x_e, epsilon_small);
    % x_h = max(x_h, epsilon_small);
    % 
    % % Compute the first term (independent of k*T) 
    % term_e = h^2 * params.dQW * nA / (VA * (4 * pi * params.m0 * params.meA));
    % term_h = h^2 * params.dQW * nA / (VA * (4 * pi * params.m0 * params.mhhA));
    % 
    % % Compute EfA
    % EfA = Eg_A + params.etaA * ( ...
    %     term_e + k * T * log(1 - exp(-x_e)) + ...
    %     term_h + k * T * log(1 - exp(-x_h)) ...
    % );

    EfB = calcEfB(T,RT,nB,params);

    % Absorption coefficients
    alpha_tm = calcalfaatm(T, RT, params.alfa_atmRT, params.alfa_atmRTdT);
    alpha_abs = calcalfaabs(T, RT, params.alfa_absRT, params.alfa_absRTdT, params.alfa_absRTdT2);
    alpha_bm = calcalfaabm(T, RT, params.alfa_abmRT, params.alfa_abmRTdT);

    eta_inj = params.eta_injxy;
    rth = calcrth(T, RT, params.rthRT, params.drthRTdT);
    tau_th = params.tau_th;

    % Drive current and internal current
    ISCH = current_func(t);  % Interpolated current
    %ISCH = current_func;  % Interpolated current
    i_inj = eta_inj * ISCH / q;

    % Barrier dynamics
    i_cap = nB / tau_cap;
    i_sp_B = nB / tau_sp_B;

    % QW dynamics
    i_esc = params.A0A * VA * T^2 / q / params.dA * exp(-Eg_B / params.eta_esc / k / T) * (exp(EfA / params.eta_esc / k / T) - 1);
    i_leak = params.A0B * VB * T^2 / q / params.dB * exp(-Eg_C / params.eta_leak / k / T) * (exp(EfB / params.eta_leak / k / T) - 1);
    if isnan(i_esc), i_esc = 0; end
    if isnan(i_leak), i_leak = 0; end

    i_sp_A = A * nA + B * nA^2 / VA + C * nA^3 / VA^2;

    % Gain
    g = (a0 * log((nA / VA + N_S) / (N0 + N_S)) - a1 * (lambda - lambdag) - a2 * (lambda - lambdag)^2) * ...
        (1 / (1 + epsilon * S * params.Gamma / VA));
    g = real(g);
    i_st = params.Gamma * params.v_g * g * S;

    % Photon dynamics
    i_abs = alpha_abs * S;
    i_tm = alpha_tm * S;
    isp = params.beta * B * nA^2 / VA;
    i_bm = alpha_bm * S;

    USCH = (nA * EfA + nB * EfB) / q / (nA + nB);

    rm = calcrm(T, RT, params.rm, params.drm, params.drm2, params.drm3);
    U = USCH*q;
    U = USCH*q+i_inj*rm*q;
    %U = USCH*q/eta_inj;
    %imped = USCH*q/ISCH

    % Temperature dynamics
    g_gen = ISCH * U - i_tm * h * c / lambda;
    g_diss = (T - params.TAMB) / rth;
    dTdt = (g_gen - g_diss) / tau_th * rth;
    if isnan(dTdt),dTdt=0;end
    %dTdt = 0;

    % Rate equations
    dnBdt = i_inj - i_sp_B - i_cap + i_esc - i_leak;
    dnAdt = i_cap - i_sp_A - i_esc - i_st;
    dSdt = isp + i_st - i_abs - i_bm - i_tm;

    % Combine derivatives
    dydt = [dnBdt; dnAdt; dSdt; dTdt];

    if ~isreal(dydt)
        aaaa = 0;
    end

    if i_tm>3.3022e+16
        dfdsfsdfds = 0;
    end

end

% Auxiliary functions

function EfB = calcEfB(T,RT,nB,params)
 
    Eg_B = calcEgB(T, RT, params.EgBRT, params.dEgBRTdT);

    VB = params.dB * params.Aaperture;  % Volume of barrier region
    % Effective density of states
    NCT = 2 * (2 * pi * params.m0 * params.meB * params.k * T / params.h^2).^(3/2);
    NVT = 2 * (2 * pi * params.m0 * params.mhhB * params.k * T / params.h^2).^(3/2);

    %Liczenie EfB
    var1C = nB ./ VB ./ NCT;
    var1V = nB ./ VB ./ NVT;
    var2C = (3 * sqrt(pi) * nB ./ VB ./ (4 * NCT)).^(2/3);
    var2V = (3 * sqrt(pi) * nB ./ VB ./ (4 * NVT)).^(2/3);

    el1C = log(var1C) ./ (1 - var1C.^2);
    el2C = var2C ./ (1 + (0.24 + 1.08 * var2C).^(-2));
    el1V = log(var1V) ./ (1 - var1V.^2);
    el2V = var2V ./ (1 + (0.24 + 1.08 * var2V).^(-2));

    EfB = Eg_B + params.etaB .* params.k .* T .* (el1C + el2C + el1V + el2V);


end

function EfA = calcEfA(T,RT,nA,params)
    Eg_A = calcEgA(T, RT, params.EgART, params.dEgARTdT);
    VA = params.dA * params.Aaperture;  % Volume of active region (QWs)
    VB = params.dB * params.Aaperture;  % Volume of barrier region
    EfA = Eg_A + params.etaA .*  ( ...
        params.k .* T .*log(exp(params.h^2 * params.dQW .* nA ./ VA ./ (4 * pi * params.k * T * params.m0 * params.meA)) - 1) + ...
        params.k .* T .*log(exp(params.h^2 * params.dQW .* nA ./ VA ./ (4 * pi * params.k * T * params.m0 * params.mhhA)) - 1) ...
    );
end

function rm = calcrm(T, RT, rmRT, drm, drm2, drm3)
    rm = rmRT + (T - RT) * drm + 0.5 * (T - RT).^2 * drm2 + (1 / 6) * (T - RT).^3 * drm3;
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
    lambdag = lambdagRT + (T - RT) * dlambdagdT + nA / VA * dlambdagdnA + 0.5 * (nA / VA)^2 * dlambdagdnA2;
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

function U = calcUdiode(current_signal,T,RT,nA,nB,params)
current_signal = current_signal(:);

EfA = calcEfA(T,RT,nA,params);
EfB = calcEfB(T,RT,nB,params);
USCH = (nA .* EfA + nB .* EfB) ./ params.q ./ (nA + nB);
i_inj = params.eta_injxy * current_signal' / params.q;
rm = calcrm(T, RT, params.rm, params.drm, params.drm2, params.drm3);
U = USCH.*params.q+i_inj.*rm.*params.q;

end

% Define the exp_safe function
function y = exp_safe(x)
    % Threshold values based on machine precision
    max_log = log(realmax);
    min_log = log(realmin);

    % Clip x to prevent numerical issues
    x_clipped = min(max(x, min_log), max_log);

    % Compute exponential
    y = exp(x_clipped);

    % Handle cases where x > max_log
    y(x > max_log) = realmax;

    % Handle cases where x < min_log
    y(x < min_log) = 0;
end