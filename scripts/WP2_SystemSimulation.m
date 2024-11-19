%% WP2_SystemSimulation.m
% This script demonstrates a complete FBG sensor system simulation
% including code generation, FBG reflection, signal processing, and detection.

% Initialize workspace
clear all;
close all;
clc;

% Ensure all folders are in the path
AddAllSubfolders();

%% Configure Simulation Parameters
% System parameters
scenariusz = 'kody1';    % Measurement scenario: 'tof', 'kody1', 'kody2', 'kody-nosna', 'tof-nosna'
kod = 'kasami';          % Code type: 'gold', 'walsh', 'kasami', 'prbs', 'randi'
mode = 'unipolar';       % Signal mode: 'unipolar' or 'bipolar'
p = 8;                   % Code length power (2^p-1)
Nb = 8;                  % Samples per bit
Fsample = 4e9;           % Sampling frequency [Hz]
K = 3;                   % Number of averages
MM_samples = 0;          % Moving mean window size

% SNR and noise parameters
noise = struct();
noise.type = 'nep';      % Noise type: 'awgn-snr', 'snr-relative', 'nep', 'true'
noise.SNR = 10;          % SNR for AWGN noise
noise.NEP = 15e-12;      % Noise equivalent power [W/√Hz]
noise.Fn = 1;            % Noise figure

% Laser parameters
laser = struct();
laser.lasers_mode = 'number';      % Wavelength generation mode: 'pm' or 'number'
laser.step = 16;                   % Number of wavelengths
laser.wavelength_range = [1534.5, 1536.5]; % Wavelength range [nm]
laser.FWHM = 2.4;                  % Full width at half maximum [pm]
laser.power = 10^(-12/10)*1e-3;    % Optical power [W]
laser.shape = 'lorentz';           % Spectral shape: 'lorentz' or 'gauss'

% Photodetector parameters
pd = struct();
pd.A = 0.8;              % Responsivity [A/W]
pd.BW = 2e9;             % Bandwidth [Hz]
pd.gain = 1;             % Linear gain
pd.RL = 50;              % Load resistance [Ω]
pd.Idark = 0.05e-12;     % Dark current [A]

% Fiber parameters
fiber = struct();
fiber.alpha =  ; % Fiber attenuation [linear]

% FBG parameters
fbg = struct();
fbg.N_s = 5;             % Number of gratings
D_spacing = 20;          % Spacing between gratings [m]
D_start = 400;           % Position of first grating [m]
fbg.D = D_start:D_spacing:(D_start+D_spacing*(fbg.N_s-1)); % Grating positions [m]
fbg.periods = 530.7300e-9 * ones(1, fbg.N_s);              % Grating periods [m]
fbg.deltaneffs = 0.3e-4 * ones(1, fbg.N_s);                % Index modulation depth
fbg.grating_lengths = 2.5e-3 * ones(1, fbg.N_s);           % Grating length [m]

% Modulation parameters (for carrier-based scenarios)
modulation = struct();
modulation.method = 'amssb';       % Modulation method
modulation.A0 = 0.1;               % Amplitude for bit 0
modulation.Fn_s = (20:90:(20+90*(laser.step-1)))*1e6; % Carrier frequencies [Hz]

% Signal processing parameters
initial_filter = struct();
initial_filter.type = 'none';      % Initial filter type: 'rrcos', 'matched', 'none'
initial_filter.alpha = 0;          % Roll-off factor for RRCOS

denoise_m = struct();
denoise_m.method = 'sgolay';       % Denoising method: 'wavelet', 'tvd', 'sgolay', 'none'
denoise_m.lambda = 0.3;            % Regularization parameter for TVD
denoise_m.order = 3;               % Order for Savitzky-Golay filter
denoise_m.framelen = 11;           % Frame length for Savitzky-Golay filter
denoise_m.iterations = 100;        % Number of iterations for TVD

%% Run Basic Simulation
fprintf('Running simulation with %d FBGs, %s coding\n', ...
    fbg.N_s, kod);

% Run system simulation
results = simulationFunction(scenariusz, kod, mode, p, Nb, Fsample, laser, ...
    K, fbg, MM_samples, modulation, initial_filter, denoise_m, pd, fiber, noise);

%% Analyze Results
% Extract Bragg wavelengths
R_org = cell2mat(results.R_simulated');
R_rec = results.R_received;

% Calculate original and recovered Bragg wavelengths
lB_org = zeros(1, fbg.N_s);
lB_rec = zeros(1, fbg.N_s);

for i = 1:fbg.N_s
    % Original Bragg wavelength
    lB_org(i) = 2 * 1.447 * fbg.periods(i) * (1 + fbg.deltaneffs(i)/1.447) * 1e9; % [nm]
    
    % Recovered Bragg wavelength (using alpha = 0.3 for weighted average)
    y_temp = R_rec(i, :);
    mx = max(y_temp) * 0.3;
    I = y_temp >= mx;
    y_temp = y_temp(I);
    x_temp = results.lambdas(I) * 1e9; % [nm]
    lB_rec(i) = sum(y_temp .* x_temp) / sum(y_temp);
end

% Calculate wavelength detection errors
lB_errors = lB_rec - lB_org;
mean_error = mean(abs(lB_errors));
std_error = std(lB_errors);

% Display results
fprintf('\n=== Wavelength Detection Results ===\n');
fprintf('Grating #\tOriginal λᵦ [nm]\tRecovered λᵦ [nm]\tError [pm]\n');
fprintf('---------------------------------------------------------\n');
for i = 1:fbg.N_s
    fprintf('%d\t\t%.4f\t\t%.4f\t\t%.2f\n', ...
        i, lB_org(i), lB_rec(i), lB_errors(i)*1000);
end
fprintf('---------------------------------------------------------\n');
fprintf('Mean Absolute Error: %.2f pm\n', mean_error*1000);
fprintf('Standard Deviation: %.2f pm\n', std_error*1000);

%% Test Successive Interference Cancellation
% Apply SIC to improve detection in case of overlapping gratings
fprintf('\n=== Testing Successive Interference Cancellation ===\n');

% Get processed data
xcv_data = results.dane_po_przetworzeniu_without_envelope;

% Apply SIC to first channel
[data_sic, detected] = successiveInterferenceCancellation(xcv_data(1,:), ...
    results.dane_upsamplowane, fbg.N_s);

% Process results with and without SIC
fprintf('\nDetection positions without SIC:\n');
[~, R_no_sic] = findGratings(xcv_data, results.dane_upsamplowane, ...
    results.length_domain, fbg.N_s, results.D_s);

fprintf('\nDetection positions with SIC:\n');
[~, R_sic] = findGratings(data_sic, results.dane_upsamplowane, ...
    results.length_domain, fbg.N_s, results.D_s);

% Calculate improvement
improvement = 100 * (sum(R_sic(:)) - sum(R_no_sic(:))) / sum(R_no_sic(:));
fprintf('\nSIC Improvement: %.2f%%\n', improvement);

%% Plot Comparative Results
% Plot original vs. SIC-processed signals
figure('Name', 'SIC Comparison', 'Position', [100, 100, 1000, 600]);

subplot(2, 1, 1);
plot(results.length_domain, xcv_data(1, :), 'LineWidth', 1.5);
hold on;
plot(results.length_domain(results.D_s), xcv_data(1, results.D_s), 'ro', 'MarkerSize', 8);
title('Signal Before SIC');
xlabel('Position [m]');
ylabel('Amplitude');
xlim([min(results.D_s)-50, max(results.D_s)+50]);
grid on;

subplot(2, 1, 2);
plot(results.length_domain, data_sic, 'LineWidth', 1.5);
hold on;
plot(results.length_domain(results.D_s), data_sic(results.D_s), 'ro', 'MarkerSize', 8);
title('Signal After SIC');
xlabel('Position [m]');
ylabel('Amplitude');
xlim([min(results.D_s)-50, max(results.D_s)+50]);
grid on;

%% Test Different Denoising Methods
% Compare various denoising methods
methods = {'none', 'wavelet', 'tvd', 'sgolay'};
results_denoising = cell(1, length(methods));

figure('Name', 'Denoising Methods Comparison', 'Position', [100, 100, 1200, 800]);

for i = 1:length(methods)
    % Set denoising method
    denoise_m.method = methods{i};
    
    % Apply denoising
    denoised_data = denoiseSignal(results.dane_w_kanale, denoise_m);
    
    % Process denoised data
    processed_data = calcXcov(denoised_data, results.dane_upsamplowane, true);
    
    % Remove envelope
    processed_without_envelope = processed_data - ...
        movmean(processed_data', results.D_s(2) - results.D_s(1))';
    
    % Calculate detection quality
    [~, R_denoised] = findGratings(processed_without_envelope, ...
        results.dane_upsamplowane, results.length_domain, fbg.N_s, results.D_s);
    
    % Store results
    results_denoising{i} = R_denoised;
    
    % Plot results
    subplot(length(methods), 1, i);
    plot(results.length_domain, processed_without_envelope(1, :), 'LineWidth', 1.5);
    hold on;
    plot(results.length_domain(results.D_s), processed_without_envelope(1, results.D_s), ...
        'ro', 'MarkerSize', 8);
    title(['Denoising Method: ' upper(methods{i})]);
    xlabel('Position [m]');
    ylabel('Amplitude');
    xlim([min(results.D_s)-50, max(results.D_s)+50]);
    grid on;
end

% Compare denoising methods performance
psnr_values = zeros(1, length(methods));
for i = 1:length(methods)
    psnr_values(i) = mean(sum(results_denoising{i}, 2));
end

figure('Name', 'Denoising Methods Performance', 'Position', [100, 100, 800, 400]);
bar(categorical(methods), psnr_values);
title('Denoising Methods Performance Comparison');
xlabel('Method');
ylabel('Average Detection Amplitude');
grid on;

%% Save Results
% Save key results for later analysis
save('WP2_simulation_results.mat', 'results', 'lB_org', 'lB_rec', ...
    'lB_errors', 'results_denoising', 'psnr_values');

fprintf('\nSimulation completed. Results saved to WP2_simulation_results.mat\n');
