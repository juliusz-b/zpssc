function [out] = filterCarrierData(data_in, data_original, Fsample, Nb, L, modulation)

%%% Projektowanie filtrow dla poszczegolnych sygnalow
out = [];
for i=1:length(modulation.Fn_s)
   fn = modulation.Fn_s(i);
   
%    B = (2*Fsample/Nb)*1e-6;B=ceil(B);B=B*1e6;
%    
%    designed_filter = designfilt('bandpassiir','FilterOrder',30, ...
%          'HalfPowerFrequency1',fn-B/2,'HalfPowerFrequency2',fn+B/2, ...
%          'SampleRate',Fsample);
%    designed_filter = designfilt('bandpassfir','FilterOrder',50, ...
%          'CutoffFrequency1',fn-B/2,'CutoffFrequency2',fn+B/2, ...
%          'SampleRate',Fsample);
%     %fvtool(designed_filter)
%     temp = filter(designed_filter,data_in);
%     %temp = bandpass(data_in,[fn-B/2 fn+B/2],Fsample);
%     plotSpectra([data_in;temp], Fsample);
%     
%     carrier = cos(2*pi*fn*(0:(Nb)-1)/Fsample);
%     
%     carrier_sig = cos(2*pi*fn*(0:(Nb*L)-1)/Fsample);
%     
%     xd = conv(temp, carrier);
%     xddd = cos(2*pi*fn*(0:(length(temp))-1)/Fsample);
%     
%     sig_mixed = temp.*xddd;
%    
%    designed_filter = designfilt('lowpassiir','FilterOrder',30, ...
%          'PassbandFrequency',B/2,'PassbandRipple',0.9, ...
%          'SampleRate',Fsample);
%     %fvtool(designed_filter)
%     lool = filter(designed_filter,sig_mixed);
%     lool = lowpass(sig_mixed,B/2-1e6,Fsample);
    temp = demod(data_in,fn,Fsample, modulation.method);
    
%     D = finddelay(data_original(i,:), lool);
%     
%     figure();plot(data_original(i,:));hold on;yyaxis right;plot(lool(D+1:end))
%     if i==2
%         plotSpectra(temp, Fsample);
%         figure();plot(temp);
%         figure();plot(xcov(temp, data_original(i,:)));
%     end
    
    out = [out; temp];
    
    
%     figure();plot(temp./carrier_sig);
%     plotSpectra(temp./carrier_sig, Fsample);
   
end

end

