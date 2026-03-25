function [data_out] = addNoise(data_in, laser, noise, pd)

%%DATA IN IS IN CURRENT DOMAIN (ELECTRICAL!)

data_out = data_in;

%snr is given in electrical domain, so:
T = 25;
T = T + 273.15;
e = 1.602176634*1e-19;
if ~isreal(data_in) && ~strcmp(noise.type,'awgn-snr')
    error('Input data is complex. Cannot add noise - switch to AWGN mode and verify input variables.')
end

snrs = [];
for i=1:size(data_in, 1)
    switch noise.type
        case 'awgn-snr'
            data_out(i,:) = awgn(data_in(i,:), noise.SNR, 'measured');
        case 'snr-relative'
            I2_rec = laser.power^2/2*pd.A^2; % laser power multiplied by responsivity
            N2 = I2_rec./(10.^(noise.SNR/10));
            N0 = sqrt(N2);
            data_out(i,:) = data_in(i,:) + N0*randn(size(data_in(i,:)));
        case 'nep'
            I2_rec = mean(data_in(i,:).^2);

            %% Electrical domain input
            P_ele = mean(data_in(i,:).^2);
            i_rms = noise.NEP*sqrt(pd.BW)*pd.A;
            PN_ele = (i_rms)^2;
            snr_ele = 10*log10(P_ele/PN_ele);
            P_opt = sqrt(mean((data_in(i,:)/pd.A).^2));
            PN_opt = noise.NEP*sqrt(pd.BW);
            snr_opt = 10*log10(P_opt/PN_opt);

            N2 = noise.NEP^2*pd.BW*pd.A^2;
            N2_optical = noise.NEP*sqrt(pd.BW);
            N2_rms_electrical = N2_optical*pd.A;N2_sqrd_electrical = N2_rms_electrical^2;
            N0 = sqrt(N2);
            data_out(i,:) = data_in(i,:) + sqrt(PN_ele)*randn(size(data_in(i,:)));
            snrs(i) = snr_ele;
        case 'true'
            I2_thermal = (4*physconst('Boltzmann')*T/pd.RL*noise.Fn*pd.BW);
            I_dark = 0.05 * 1e-12;
            I2_electrical = mean(data_in(i,:).^2);
            N2Total = (2*e*pd.BW*(I_dark+sqrt(I2_electrical))+I2_thermal);
            N0Total = sqrt(N2Total);
            data_out(i,:) = data_in(i,:) + N0Total*randn(size(data_in(i,:)));
            snrs(i) = 10*log10(I2_electrical/N2Total);
        otherwise
            error('Unknown noise type');
    end
end

if ~evalin('base', 'exist(''SILENT_MODE'',''var'') && SILENT_MODE')
    disp(minmax(snrs))
end
end

