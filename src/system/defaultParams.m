function params = defaultParams()
% DEFAULTPARAMS Return default simulation parameters
%
% params = defaultParams()
%
% Returns a struct with all parameters set to sensible defaults.
% Modify individual fields before passing to runSimulation().

params.scenariusz = 'kody2';
params.kod = 'kasami';
params.mode = 'unipolar';
params.p = 8;
params.Nb = 8;
params.Fsample = 4e9;
params.K = 1;
params.MM_samples = 0;
params.show_plots = false;  % Set to true to display plots after simulation

% Noise
params.noise.type = 'nep';
params.noise.SNR = 10;
params.noise.NEP = 15e-12;
params.noise.Fn = 1;

% Laser
params.laser.lasers_mode = 'number';
params.laser.step = 16;
params.laser.wavelength_range = [1534.5, 1536.5];
params.laser.FWHM = 2.4;
params.laser.power = 10^(-12/10)*1e-3;
params.laser.shape = 'lorentz';

% Photodetector
params.pd.A = 0.8;
params.pd.BW = 250e6;
params.pd.gain = 1;
params.pd.RL = 50;
params.pd.Idark = 0.05e-12;

% Fiber
params.fiber.alpha = 1;

% FBG
params.fbg.N_s = 5;
D_spacing = 20;
D_start = 400;
params.fbg.D = D_start:D_spacing:(D_start + D_spacing*(params.fbg.N_s-1));
params.fbg.periods = 530.7300e-9 * ones(1, params.fbg.N_s);
params.fbg.deltaneffs = 0.3e-4 * ones(1, params.fbg.N_s);
params.fbg.grating_lengths = 3e-3 * ones(1, params.fbg.N_s);
params.fbg.model = 'tanh';
params.fbg.max_bounces = 1;
params.fbg.lambda_shifts = zeros(1, params.fbg.N_s);

% Modulation
params.modulation.method = 'amssb';
params.modulation.A0 = 0.1;
params.modulation.Fn_s = (20:90:(20+90*(params.laser.step-1)))*1e6;

% Initial filter
params.initial_filter.type = 'none';
params.initial_filter.alpha = 0;

% Denoising
params.denoise_m.method = 'none';
params.denoise_m.lambda = 0.3;
params.denoise_m.order = 3;
params.denoise_m.framelen = 11;
params.denoise_m.iterations = 100;

end
