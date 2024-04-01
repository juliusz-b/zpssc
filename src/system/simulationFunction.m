function [outputResults] = simulationFunction(scenariusz, kod, mode, p, Nb, Fsample, laser, K, SNR, fbg, MM_samples, modulation, initial_filter, denoise_m, pd, fiber, noise)
% SIMULATIONFUNCTION Main simulation function for FBG sensing system
%
% [outputResults] = simulationFunction(scenariusz, kod, mode, p, Nb, Fsample, laser, K, SNR, fbg, MM_samples, modulation, initial_filter, denoise_m, pd, fiber, noise)
%
% Inputs:
%   scenariusz - Type of scenario: 'tof', 'kody1', 'kody2', 'kody-nosna', 'tof-nosna'
%   kod - Code type: 'gold', 'walsh', 'kasami', 'prbs', 'randi'
%   mode - Signal mode: 'unipolar' or 'bipolar'
%   p - Power of sequence length (2^p-1)
%   Nb - Samples per bit
%   Fsample - Sampling frequency [Hz]
%   laser - Laser parameters structure
%   K - Number of averages
%   SNR - Signal-to-noise ratio [dB]
%   fbg - FBG parameters structure
%   MM_samples - Moving mean window size
%   modulation - Modulation parameters
%   initial_filter - Initial filter parameters
%   denoise_m - Denoising method parameters
%   pd - Photodetector parameters
%   fiber - Fiber parameters
%   noise - Noise parameters
%
% Outputs:
%   outputResults - Structure containing simulation results

% Calculate code length
L = 2^p-1;

% Generate wavelengths for laser
lambdas = generateLambdas(laser);

% Simulate laser source
[simulation_range, lasers_spectra] = generateLasers(laser);

% Generate FBG reflection spectra
[R_s, R_simulated] = generateAllGratings(lambdas, fbg, simulation_range);

% Determine number of wavelength channels based on scenario
switch scenariusz
    case 'tof'
        U = length(lambdas); % U equals 1 for impulse
        L = 1;
        Nb = 1;
    case {'kody1', 'kody2', 'kody-nosna', 'tof-nosna'}
        U = length(lambdas); % Number of wavelength channels
    otherwise
        error('Invalid scenario type!');
end

% Generate distance values and domains
[D_s, L_max, time_domain, length_domain] = generateValues(L, Fsample, Nb, fbg.D, U, fbg.N_s, false);

% Select mode function (unipolar/bipolar) and power parameter
[BipUni, power] = selectMode(mode);

% Generate spreading codes
if contains(scenariusz, 'tof')
    % For TOF, use impulse sequences
    kody = BipUni([ones(U, 1), zeros(U, L-1)]);
else
    % For other scenarios, use spreading codes
    kody = BipUni(genSpreadCodes(L, U, kod));
end

% Upsample codes to match required sample rate
dane_upsamplowane = myUpsample(kody, Nb);

% Generate initial filter
[filter_data] = generateInitialFilter(initial_filter, Fsample, Nb, 16);

% Apply initial filtering if necessary
if ~strcmp(initial_filter.type, 'matched')
    dane_wstepnie_filtrowane = filterMySignal(dane_upsamplowane, filter_data);
else
    dane_wstepnie_filtrowane = dane_upsamplowane;
end

% Apply carrier modulation for relevant scenarios
if contains(scenariusz, 'nosna')
    % Generate carrier signals
    sygnaly_nosne = generateCarrierSignals(L, U, Nb, Fsample, modulation);
    
    % Apply modulation to carriers
    dane_na_nosnych = moveSignalToCarrier(sygnaly_nosne, dane_wstepnie_filtrowane, mode, modulation, Fsample);
    
    % Add gratings to modulated signals
    dane_z_siatkami = addAllGratings(dane_na_nosnych, fbg.N_s, D_s, R_s, L_max, R_simulated, lasers_spectra, simulation_range, fbg, fiber);
else
    % Add gratings to filtered signals without modulation
    dane_z_siatkami = addAllGratings(dane_wstepnie_filtrowane, fbg.N_s, D_s, R_s, L_max, R_simulated, lasers_spectra, simulation_range, fbg, fiber);
end

% Process signals based on scenario
switch scenariusz
    case 'tof'
        % Time of flight measurement
        % Add noise and apply photodetector response
        dane_w_kanale = noiseAndMeanDataInChannel(dane_z_siatkami, SNR, K, laser, noise, pd).^power;
        dane_w_kanale = lowpass(dane_w_kanale, pd.BW, Fsample);
        dane_w_kanale = dane_w_kanale * pd.A * pd.gain;
        
        % Apply filtering and denoising
        dane_w_kanale = filterMySignal(dane_w_kanale, filter_data);
        dane_w_kanale = denoiseSignal(dane_w_kanale, denoise_m);
        
        % Output processed data
        dane_po_przetworzeniu = dane_w_kanale;
        
    case 'kody1'
        % Code division sensor system - separate channels
        % Convert optical power to electrical current
        dane_w_kanale = dane_z_siatkami * pd.A;
        
        % Add noise and apply electrical filtering
        dane_w_kanale = noiseAndMeanDataInChannel(dane_w_kanale, SNR, K, laser, noise, pd).^power;
        
        % Apply bandwidth limitation
        for i = 1:size(dane_w_kanale, 1)
            dane_w_kanale(i,:) = lowpass(dane_w_kanale(i,:), pd.BW, Fsample);
        end
        
        % Apply photodetector gain
        dane_w_kanale = dane_w_kanale * pd.gain;
        
        % Apply filtering and denoising
        dane_w_kanale = filterMySignal(dane_w_kanale, filter_data);
        dane_w_kanale = denoiseSignal(dane_w_kanale, denoise_m);
        
        % Calculate cross-covariance for detection
        dane_po_przetworzeniu = calcXcov(dane_w_kanale, dane_upsamplowane, true);
        
    case 'kody2'
        % Code division sensor system - combined channels
        % Sum optical powers from all channels
        dane_w_kanale = sum(dane_z_siatkami, 1);
        
        % Convert to electrical domain
        dane_w_kanale = dane_w_kanale * pd.A;
        
        % Add noise
        dane_w_kanale = noiseAndMeanDataInChannel(dane_w_kanale, SNR, K).^power;
        
        % Apply bandwidth limitation
        for i = 1:size(dane_w_kanale, 1)
            dane_w_kanale(i,:) = lowpass(dane_w_kanale(i,:), pd.BW, Fsample);
        end
        
        % Apply photodetector gain
        dane_w_kanale = dane_w_kanale * pd.gain;
        
        % Apply filtering and denoising
        dane_w_kanale = filterMySignal(dane_w_kanale, filter_data);
        dane_w_kanale = denoiseSignal(dane_w_kanale, denoise_m);
        
        % Calculate cross-covariance for detection
        dane_po_przetworzeniu = calcXcov(dane_w_kanale, dane_upsamplowane, true);
        
    case 'kody-nosna'
        % Code division with carrier modulation
        % Add noise to summed signals
        dane_w_kanale = noiseAndMeanDataInChannel(sum(dane_z_siatkami, 1), SNR, K).^power;
        dane_w_kanale = lowpass(dane_w_kanale, pd.BW, Fsample);
        dane_w_kanale = dane_w_kanale * pd.A * pd.gain;
        
        % Filter carrier data
        dane_po_filtrze = filterCarrierData(dane_w_kanale, dane_upsamplowane, Fsample, Nb, L, modulation);
        
        % Apply filtering and denoising
        dane_po_filtrze = filterMySignal(dane_po_filtrze, filter_data);
        dane_w_kanale = denoiseSignal(dane_po_filtrze, denoise_m);
        
        % Calculate cross-covariance for detection
        dane_po_przetworzeniu = abs(calcXcov(dane_w_kanale, dane_upsamplowane, true));
        
    case 'tof-nosna'
        % Time of flight with carrier modulation
        % Add noise to summed signals
        dane_w_kanale = noiseAndMeanDataInChannel(sum(dane_z_siatkami, 1), SNR, K).^power;
        dane_w_kanale = lowpass(dane_w_kanale, pd.BW, Fsample);
        dane_w_kanale = dane_w_kanale * pd.A * pd.gain;
        
        % Filter carrier data
        dane_po_filtrze = filterCarrierData(dane_w_kanale, dane_upsamplowane, Fsample, Nb, L, modulation);
        
        % Apply filtering and denoising
        dane_w_kanale = filterMySignal(dane_po_filtrze, filter_data);
        dane_w_kanale = denoiseSignal(dane_w_kanale, denoise_m);
        
        % For TOF, processed data is the channel data
        dane_po_przetworzeniu = dane_w_kanale;
        
    otherwise
        error('Invalid scenario type!');
end

% Apply moving mean if specified
if MM_samples > 0 && ~strcmp(scenariusz, 'tof') && Nb > 1
    dane_po_przetworzeniu = dane_po_przetworzeniu - movmean(dane_po_przetworzeniu, MM_samples);
end

% Remove correlation envelope
dane_po_przetworzeniu_without_envelope = dane_po_przetworzeniu - movmean(dane_po_przetworzeniu', D_s(2) - D_s(1))';

% Find gratings and calculate reflectance
[~, R_received] = findGratings(dane_po_przetworzeniu_without_envelope, dane_upsamplowane, length_domain, fbg.N_s, D_s);
R_received = R_received.^(1/power) / laser.power / pd.A;

% Upsample and denoise FBG reflectance signals
[new_lmb, R_filtered, lmbB] = upsampleAndDenoiseGratings(R_received, fbg, lambdas);

% Prepare output structure
outputResults = struct();
outputResults.lambdas = lambdas;
outputResults.R_received = R_received;
outputResults.R_filtered = R_filtered;
outputResults.new_lmb = new_lmb;
outputResults.lmbB = lmbB;
outputResults.dane_po_przetworzeniu = dane_po_przetworzeniu;
outputResults.dane_po_przetworzeniu_without_envelope = dane_po_przetworzeniu_without_envelope;
outputResults.dane_w_kanale = dane_w_kanale;
outputResults.dane_upsamplowane = dane_upsamplowane;
outputResults.length_domain = length_domain;
outputResults.time_domain = time_domain;
outputResults.D_s = D_s;
outputResults.R_simulated = R_simulated;
outputResults.simulation_range = simulation_range;

% Display results if not suppressed
plotGratingsSnr(lambdas, R_s, SNR, fbg);
plotGratingReflectanceSurf(lambdas, R_simulated, R_received, length_domain, ...
                          dane_po_przetworzeniu_without_envelope, dane_w_kanale, ...
                          dane_upsamplowane, fbg.N_s, D_s, R_simulated, ...
                          simulation_range, new_lmb, R_filtered);
plotXcovs(dane_po_przetworzeniu_without_envelope, length_domain);
end
