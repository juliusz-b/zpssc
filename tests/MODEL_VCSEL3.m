% Large-Signal VCSEL Model with Signal Vector Input

%REFERENCE:A. Grabowski, J. Gustavsson, Z. S. He, and A. Larsson, “Large-Signal Equivalent Circuit for Datacom VCSELs,” Journal of Lightwave Technology,
% vol. 39, no. 10. Institute of Electrical and Electronics Engineers (IEEE), pp. 3225–3236, May 15, 2021. doi: 10.1109/jlt.2021.3064465.

Ibias = 40e-3;      % Bias current (A)
Imod = 5e-3;
UF = 4;
SR = 40e9;
SRCAP = SR;
BR = SRCAP/UF;
N = 2e4;
Mary = 4;

genPAM4_Fudan;

%NLED = length(x_in);
%time_vector = 0:1/SR:(1/SR)*(NLED-1);

%time_vector = linspace(0, 100e-9, 10000);  % Time vector in seconds
%current_signal = Ibias + (Ibias / 2) * sin(2 * pi * 1e9 * time_vector);  % Example signal
current_signal = Ibias + (Imod) * x_in;  % Example signal

%VCSELparams1;

%syg = VCSELtransmission3(0*[x_in; x_in; x_in; x_in],SR,Ibias,Imod);
syg = VCSELtransmission3([x_in],SR,Ibias,Imod);

L_PROBEK = length(syg);
recPAM4_Fudan_script;


% Plot results
% figure;
% subplot(3, 1, 1);
% plot(t, P_opt*1e3);
% xlabel('Time (s)');
% ylabel('Optical Power (mW)');
% title('Optical Power');
% 
% subplot(3, 1, 2);
% plot(t, nA);
% xlabel('Time (s)');
% ylabel('Carrier Concentration in QWs (m^{-3})');
% title('Carrier Concentration in QWs');
% 
% subplot(3, 1, 3);
% plot(t, T);
% xlabel('Time (s)');
% ylabel('Temperature (K)');
% title('Internal Temperature');
