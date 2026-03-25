function [out] = runSimulation(params)
% RUNSIMULATION Run FBG CDM simulation - clean pipeline, no plots
%
% [out] = runSimulation(params)
%
% Input: params struct with fields:
%   .scenariusz  - 'tof', 'kody1', 'kody2', 'kody-nosna', 'tof-nosna'
%   .kod         - 'kasami', 'gold', 'walsh', 'prbs'
%   .mode        - 'unipolar' / 'bipolar'
%   .p           - Code length exponent (2^p-1)
%   .Nb          - Samples per bit
%   .Fsample     - Sampling frequency [Hz]
%   .K           - Number of averages
%   .MM_samples  - Moving mean window (0 = off)
%   .laser       - Laser parameters struct
%   .pd          - Photodetector parameters struct
%   .fiber       - Fiber parameters struct
%   .fbg         - FBG parameters struct
%   .noise       - Noise parameters struct
%   .modulation  - Modulation parameters struct
%   .initial_filter - Initial filter struct
%   .denoise_m   - Denoising method struct
%
% Output: out struct with all intermediate and final results

% Unpack
scenariusz = params.scenariusz;
kod = params.kod;
mode = params.mode;
p = params.p;
Nb = params.Nb;
Fsample = params.Fsample;
K = params.K;
MM_samples = params.MM_samples;
laser = params.laser;
pd = params.pd;
fiber = params.fiber;
fbg = params.fbg;
noise = params.noise;
modulation = params.modulation;
initial_filter = params.initial_filter;
denoise_m = params.denoise_m;

neff = 1.447;

% Code length
L = 2^p - 1;

% Wavelengths
lambdas = generateLambdas(laser);

% Laser spectra
[simulation_range, lasers_spectra] = generateLasers(laser);

% FBG reflectance spectra
[R_s, R_simulated] = generateAllGratings(lambdas, fbg, simulation_range);

% Scenario channels
switch scenariusz
    case 'tof'
        U = length(lambdas);
        L = 1; Nb = 1;
    case {'kody1', 'kody2', 'kody-nosna', 'tof-nosna'}
        U = length(lambdas);
    otherwise
        error('Unknown scenario: %s', scenariusz);
end

% Spatial/temporal domain
if isfield(fbg, 'max_bounces'); mb = fbg.max_bounces; else; mb = 1; end
[D_s, L_max, time_domain, length_domain] = generateValues(L, Fsample, Nb, fbg.D, U, fbg.N_s, false, mb);

% Mode
[BipUni, power] = selectMode(mode);

% Spreading codes
if contains(scenariusz, 'tof')
    kody = BipUni([ones(U, 1), zeros(U, L-1)]);
else
    kody = BipUni(genSpreadCodes(L, U, kod));
end

% Upsample
dane_upsamplowane = myUpsample(kody, Nb);

% Initial filter
[filter_data] = generateInitialFilter(initial_filter, Fsample, Nb, 16);
if ~strcmp(initial_filter.type, 'matched')
    dane_wstepnie_filtrowane = filterMySignal(dane_upsamplowane, filter_data);
else
    dane_wstepnie_filtrowane = dane_upsamplowane;
end

% Carrier modulation (optional)
if contains(scenariusz, 'nosna')
    sygnaly_nosne = generateCarrierSignals(L, U, Nb, Fsample, modulation);
    dane_do_siatek = moveSignalToCarrier(sygnaly_nosne, dane_wstepnie_filtrowane, mode, modulation, Fsample);
else
    dane_do_siatek = dane_wstepnie_filtrowane;
end

% FBG reflections
dane_z_siatkami = addAllGratings(dane_do_siatek, fbg.N_s, D_s, R_s, L_max, ...
    R_simulated, lasers_spectra, simulation_range, fbg, fiber);

% Design FIR lowpass for photodetector bandwidth (faster than lowpass())
b_lpf = fir1(64, pd.BW / (Fsample/2));
applyLPF = @(x) filter(b_lpf, 1, x);

% Channel processing
switch scenariusz
    case 'tof'
        dane_w_kanale = noiseAndMeanDataInChannel(dane_z_siatkami, K, laser, noise, pd).^power;
        for i = 1:size(dane_w_kanale, 1); dane_w_kanale(i,:) = applyLPF(dane_w_kanale(i,:)); end
        dane_w_kanale = dane_w_kanale * pd.A * pd.gain;
        dane_w_kanale = filterMySignal(dane_w_kanale, filter_data);
        dane_w_kanale = denoiseSignal(dane_w_kanale, denoise_m);
        dane_po_przetworzeniu = dane_w_kanale;

    case 'kody1'
        dane_w_kanale = dane_z_siatkami * pd.A;
        dane_w_kanale = noiseAndMeanDataInChannel(dane_w_kanale, K, laser, noise, pd).^power;
        for i = 1:size(dane_w_kanale, 1); dane_w_kanale(i,:) = applyLPF(dane_w_kanale(i,:)); end
        dane_w_kanale = dane_w_kanale * pd.gain;
        dane_w_kanale = filterMySignal(dane_w_kanale, filter_data);
        dane_w_kanale = denoiseSignal(dane_w_kanale, denoise_m);
        dane_po_przetworzeniu = calcXcov(dane_w_kanale, dane_upsamplowane, true);

    case 'kody2'
        dane_w_kanale = sum(dane_z_siatkami, 1) * pd.A;
        dane_w_kanale = noiseAndMeanDataInChannel(dane_w_kanale, K, laser, noise, pd).^power;
        for i = 1:size(dane_w_kanale, 1); dane_w_kanale(i,:) = applyLPF(dane_w_kanale(i,:)); end
        dane_w_kanale = dane_w_kanale * pd.gain;
        dane_w_kanale = filterMySignal(dane_w_kanale, filter_data);
        dane_w_kanale = denoiseSignal(dane_w_kanale, denoise_m);
        dane_po_przetworzeniu = calcXcov(dane_w_kanale, dane_upsamplowane, true);

    case 'kody-nosna'
        dane_w_kanale = noiseAndMeanDataInChannel(sum(dane_z_siatkami, 1), K).^power;
        for i = 1:size(dane_w_kanale, 1); dane_w_kanale(i,:) = applyLPF(dane_w_kanale(i,:)); end
        dane_w_kanale = dane_w_kanale * pd.A * pd.gain;
        dane_po_filtrze = filterCarrierData(dane_w_kanale, dane_upsamplowane, Fsample, Nb, L, modulation);
        dane_po_filtrze = filterMySignal(dane_po_filtrze, filter_data);
        dane_w_kanale = denoiseSignal(dane_po_filtrze, denoise_m);
        dane_po_przetworzeniu = abs(calcXcov(dane_w_kanale, dane_upsamplowane, true));

    case 'tof-nosna'
        dane_w_kanale = noiseAndMeanDataInChannel(sum(dane_z_siatkami, 1), K).^power;
        for i = 1:size(dane_w_kanale, 1); dane_w_kanale(i,:) = applyLPF(dane_w_kanale(i,:)); end
        dane_w_kanale = dane_w_kanale * pd.A * pd.gain;
        dane_po_filtrze = filterCarrierData(dane_w_kanale, dane_upsamplowane, Fsample, Nb, L, modulation);
        dane_w_kanale = filterMySignal(dane_po_filtrze, filter_data);
        dane_w_kanale = denoiseSignal(dane_w_kanale, denoise_m);
        dane_po_przetworzeniu = dane_w_kanale;
end

% Moving mean
if MM_samples > 0 && ~strcmp(scenariusz, 'tof') && Nb > 1
    dane_po_przetworzeniu = dane_po_przetworzeniu - movmean(dane_po_przetworzeniu, MM_samples);
end

% Remove correlation envelope
dane_bez_obwiedni = dane_po_przetworzeniu - movmean(dane_po_przetworzeniu', D_s(2) - D_s(1))';

% Grating detection
[~, R_received] = findGratings(dane_bez_obwiedni, dane_upsamplowane, length_domain, fbg.N_s, D_s);
R_received = R_received.^(1/power) / laser.power / pd.A;

% Upsample and denoise FBG reflectance
[new_lmb, R_filtered, lmbB] = upsampleAndDenoiseGratings(R_received, fbg, lambdas);

% Calculate Bragg wavelength errors
lB_org = zeros(1, fbg.N_s);
lB_rec = zeros(1, fbg.N_s);
for i = 1:fbg.N_s
    % Original lambda_B (accounting for shifts)
    period_i = fbg.periods(i);
    if isfield(fbg, 'lambda_shifts') && fbg.lambda_shifts(i) ~= 0
        period_i = period_i * (1 + fbg.lambda_shifts(i) / (2*neff*period_i));
    end
    lB_org(i) = 2 * neff * period_i * (1 + fbg.deltaneffs(i)/neff) * 1e9;

    % Recovered lambda_B via centroid
    y = R_received(i, :);
    mx = max(y) * 0.3;
    I = y >= mx;
    if any(I) && sum(y(I)) > 0
        lB_rec(i) = sum(y(I) .* lambdas(I) * 1e9) / sum(y(I));
    else
        lB_rec(i) = NaN;
    end
end

% Output struct
out.lambdas = lambdas;
out.simulation_range = simulation_range;
out.R_s = R_s;
out.R_simulated = R_simulated;
out.R_received = R_received;
out.R_filtered = R_filtered;
out.new_lmb = new_lmb;
out.lmbB = lmbB;
out.dane_po_przetworzeniu = dane_po_przetworzeniu;
out.dane_bez_obwiedni = dane_bez_obwiedni;
out.dane_w_kanale = dane_w_kanale;
out.dane_upsamplowane = dane_upsamplowane;
out.dane_z_siatkami = dane_z_siatkami;
out.length_domain = length_domain;
out.time_domain = time_domain;
out.D_s = D_s;
out.L_max = L_max;
out.lB_org = lB_org;
out.lB_rec = lB_rec;
out.lB_errors = lB_rec - lB_org;
out.MAE = mean(abs(out.lB_errors), 'omitnan') * 1000; % [pm]
out.MAX_err = max(abs(out.lB_errors)) * 1000; % [pm]

% Plots (optional)
if isfield(params, 'show_plots') && params.show_plots
    % Korelacja krzyżowa
    figure('Position', [50 50 1200 450], 'theme', 'light', 'Name', 'Korelacja krzyżowa');
    n_show = min(4, size(dane_bez_obwiedni, 1));
    colors = lines(n_show);
    for ch = 1:n_show
        plot(length_domain, dane_bez_obwiedni(ch,:), 'LineWidth', 1.2, 'Color', colors(ch,:));
        hold on;
    end
    for j = 1:fbg.N_s
        xline(length_domain(D_s(j)), '--k', sprintf('FBG%d', j), 'LineWidth', 1.2, ...
            'LabelOrientation', 'horizontal', 'FontSize', 10);
    end
    xlabel('Pozycja [m]'); ylabel('Amplituda');
    title(sprintf('Korelacja - %s, %s, p=%d, MAE=%.1f pm', scenariusz, kod, p, out.MAE));
    xlim([min(fbg.D)-30, max(fbg.D)+30]);
    grid on; set(gca, 'FontSize', 12);

    % Reflektancja + błędy
    figure('Position', [50 500 1200 450], 'theme', 'light', 'Name', 'Reflektancja i błędy');
    subplot(1,2,1);
    sim_nm = simulation_range * 1e9;
    lmb_nm = lambdas * 1e9;
    colors_fbg = lines(fbg.N_s);
    for i = 1:fbg.N_s
        plot(sim_nm, R_simulated{i}, '-', 'Color', colors_fbg(i,:), 'LineWidth', 1.5); hold on;
        plot(lmb_nm, R_received(i,:), 'o', 'Color', colors_fbg(i,:), 'MarkerSize', 5, 'MarkerFaceColor', colors_fbg(i,:));
    end
    xlabel('\lambda [nm]'); ylabel('Reflektancja');
    title('Widmo: oryginalne vs odtworzone');
    grid on; set(gca, 'FontSize', 11);

    subplot(1,2,2);
    bar(out.lB_errors * 1000, 'FaceColor', [0.85 0.33 0.1]);
    xlabel('Siatka'); ylabel('\Delta\lambda_B [pm]');
    title(sprintf('Błędy (MAE=%.1f pm)', out.MAE));
    grid on; set(gca, 'FontSize', 11);
    xticklabels(arrayfun(@(x) sprintf('FBG%d', x), 1:fbg.N_s, 'UniformOutput', false));
end

end
