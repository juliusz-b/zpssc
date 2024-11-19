function [syg] = VCSELtransmission3(x_in,SR,Ibias,Imod)
% End-to-end optimization of optical communication systems based on directly modulated lasers

%REFERENCE:Sergio Hernandez, Christophe Peucheret, Francesco Da Ros, and Darko Zibar, 
% "End-to-end optimization of optical communication systems based on directly modulated lasers," J. Opt. Commun. Netw. 16, D29-D43 (2024)

%czesc parametrow z:https://doi.org/10.1038/s41566-020-00700-y


NLED = length(x_in);
time_vector = 0:1/SR:(1/SR)*(NLED-1);

current_signal = Ibias + (Imod) * x_in;  % Example signal




params.Gamma = 0.24;
    params.g0 = 3500; % z https://doi.org/10.1038/s41566-020-00700-y
params.N0 = 2e24;
 params.epsilon = 2e23;
params.tauP = 2.6e-12;
 params.beta = 1e-3;
 params.tauN = 3.17e-9; %w papierze jest podany carrier lifetime...
params.q = 1.6e-19;
 params.V = 3.6e-17; 
params.alpha = 0.1; %wymyslilem, bo nie podali....  
params.eta0 = 0.2; %
params.h = 6.626e-34; % Planck constant (J·s
params.ng = 4; % eff refr indx
params.c = 3e8/params.ng;           % Speed of light (m/s)
params.lambda = 1550e-9;   % Wavelength (m)
params.nu = params.c/params.lambda;


% Initial Conditions

%initial_conditions = [nB0, nA0, 1e-9, RT];  % [nB, nA, S, T]
%initial_conditions = [nB0 nA0 S0 RT];  % [nB, nA, S, T]
initial_conditions = [0 0 0];  % [nB, nA, S, T]

% Interpolated current function
%current_func = @(t) interp1q(time_vector, current_signal, t);
current_func = @(t) interp1(time_vector, current_signal, t, 'linear', 'extrap');


tspan = [time_vector(1), time_vector(end)];
tic
[t, y] = ode89(@(t, y) vcsel_odes(t, y, params, current_func), tspan, initial_conditions);
toc

S = y(:, 1);
N = y(:, 2);  
Phi = y(:, 3);  




% opt out

P_opt = S * params.V * params.eta0*params.h*params.nu/(2*params.Gamma*params.tauP);

P_opt = interp1(t,P_opt,time_vector);

syg = P_opt';


% Plot results
figure;
subplot(3, 2, 1);
plot(time_vector, P_opt*1e3);
xlabel('Time (s)');
ylabel('Optical Power (mW)');
title('Optical Power');







end


% Function definitions
function dydt = vcsel_odes(t, y, params, current_func)
    % Extract variables
    S = y(1); 
    N = y(2); 
    Phi = y(3);
    

    % Extract parameters
    Gamma = params.Gamma;
    g0 = params.g0;
    N0 = params.N0;
    epsilon = params.epsilon;
    tauP = params.tauP;
    beta = params.beta;
    tauN = params.tauN;
    q = params.q;
    V = params.V; 
    alpha = params.alpha;  


    dSdt = Gamma*g0*(N-N0)*1./(1+epsilon*S)*S-S/tauP+Gamma*beta*N/tauN;
    dNdt = current_func(t)/(q*V) - N/tauN - g0*(N-N0)*1./(1+epsilon*S)*S;
    dPhidt = 0*1/2*alpha*(Gamma*g0*(N-N0) - 1/tauP    );

    % Combine derivatives
    dydt = [dSdt; dNdt; dPhidt];


end