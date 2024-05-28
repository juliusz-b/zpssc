function [data_out] = addNoise(data_in, laser, noise, pd)

%%DATA IN IS IN CURRENT DOMAIN (ELECTRICAL!)

data_out = data_in;

%snr is given in electrical domain, so:
T = 25;
T = T + 273.15;
e = 1.602176634*1e-19;
if ~isreal(data_in) && ~strcmp(noise.type,'awgn-snr')
    error('Dane sa zespolone. Nie mozna dodac szumu - zmien na AWGN, ale sprawdź dobrze zmienne wprowadzane do awgn :D')
end

snrs = [];
for i=1:size(data_in, 1)
    switch noise.type
        case 'awgn-snr'
            data_out(i,:) = awgn(data_in(i,:), noise.SNR, 'measured');
        case 'snr-relative'
            I2_rec = laser.power^2/2*pd.A^2;%moc lasera przemnozona przez responsivity
            N2 = I2_rec./(10.^(noise.SNR/10));
            N0 = sqrt(N2);
            data_out(i,:) = data_in(i,:) + N0*randn(size(data_in(i,:)));
        case 'nep'
            I2_rec = mean(data_in(i,:).^2);

            %%%%%%%%%% Przypadek, gdy wejscie jest elektryczne!
            P_ele = mean(data_in(i,:).^2);
            i_rms = noise.NEP*sqrt(pd.BW)*pd.A;
            PN_ele = (i_rms)^2;
            snr_ele = 10*log10(P_ele/PN_ele);
            P_opt = sqrt(mean((data_in(i,:)/pd.A).^2));
            PN_opt = noise.NEP*sqrt(pd.BW);
            snr_opt = 10*log10(P_opt/PN_opt);
            %%%%%%%%%%
            
%             %%%%%%%%%% Przypadek, gdy wejscie jest optyczne!
%             P_opt = sqrt(mean(data_in(i,:).^2));
%             PN_opt = noise.NEP*sqrt(pd.BW);
%             snr_opt = 10*log10(P_opt/PN_opt);
%             P_ele = mean((data_in(i,:)*pd.A).^2);
%             PN_ele = (noise.NEP*sqrt(pd.BW)*pd.A)^2;
%             I2_rec_ele = mean((data_in(i,:)*pd.A).^2);
%             snr_ele = 10*log10(P_ele/PN_ele)
%             %%%%%%%%%%

            %N2 = noise.NEP*sqrt(pd.BW)*pd.A;
            N2 = noise.NEP^2*pd.BW*pd.A^2;
            %N2 = noise.NEP*sqrt(pd.BW)/pd.A;
            N2_optical = noise.NEP*sqrt(pd.BW);
            %wzorki z https://books.google.pl/books?id=aKCgh3go-sYC&pg=PA70&lpg=PA70&dq=nep+SNR&source=bl&ots=ZFe2Imcvw1&sig=ACfU3U3oO5vhO9EEhnpo8Jw3G1C4zzKydA&hl=pl&sa=X&ved=2ahUKEwi2uLm55pf6AhUOt4sKHRSlC0EQ6AF6BAgZEAM#v=onepage&q=nep%20SNR&f=false
            N2_rms_electrical = N2_optical*pd.A;N2_sqrd_electrical = N2_rms_electrical^2;
            %koniec wzorkow
            N0 = sqrt(N2);
            data_out(i,:) = data_in(i,:) + sqrt(PN_ele)*randn(size(data_in(i,:)))/1; %dzielimy na R, zeby byl to szum "optyczny"
            snrs(i) = snr_ele;
            %sprawdzenie czy OK
%             tst =  randn(size(data_in(i,:))); P_tst = mean(tst.^2);tst_n = sqrt(PN_ele)*randn(size(data_in(i,:))) + tst;e = comm.EVM;RMSEVM=step(e,tst_n',tst');EVM=RMSEVM*0.01;SNR_EVM=20*log10(1./EVM);SNR_CAL = 10*log10(P_tst/PN_ele); disp([SNR_EVM SNR_CAL])
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

% warning('on','all')
% snrs = [];
% for i=1:size(data_in, 1)
% 
%     I2_electrical = mean(data_in(i,:).^2);
%     N0 = sqrt(I2_optical./(10.^(snr/10)));
%     N2Total = (2*e*BW*(I_dark+sqrt(I2_electrical))+I2_thermal);
%     N0Total = sqrt(N2Total);
%     %data_out(i,:) = awgn(data_in(i,:), snr, 'measured');
%     if true_noise
%         data_out(i,:) = data_in(i,:) + N0Total*randn(size(data_in(i,:)));
%     else
%         data_out(i,:) = data_in(i,:) + N0*randn(size(data_in(i,:)));
%         %data_out(i,:) = awgn(data_in(i,:), snr, 'measured');
%     end
%     snrs(i) = 10*log10(I2_electrical/N2Total);
% end
disp(minmax(snrs))
warning('off','all')
end

