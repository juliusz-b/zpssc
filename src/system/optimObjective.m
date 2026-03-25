function mae = optimObjective(x, p_map, Nb_map, nlam_map, N_s, grad_nm)
% OPTIMOBJECTIVE Objective function for GA/PSO optimization
%
% mae = optimObjective(x, p_map, Nb_map, nlam_map, N_s, grad_nm)
%
% x = [deltaneff, p_idx, Nb_idx, nlam_idx, spectral_spread_nm]

deltaneff = x(1);
p = p_map(max(1, min(length(p_map), round(x(2)))));
Nb = Nb_map(max(1, min(length(Nb_map), round(x(3)))));
n_lam = nlam_map(max(1, min(length(nlam_map), round(x(4)))));
spread_nm = x(5);

L = 2^p - 1;
Fsample = 4e9;
f_interr = Fsample / (L * Nb * n_lam * 2);

% Penalty: f_interrogation < 1 kHz
if f_interr < 1000
    mae = 1e4;
    return;
end

params = defaultParams();
params.show_plots = false;
params.p = p;
params.Nb = Nb;
params.laser.step = n_lam;
params.modulation.Fn_s = (20:90:(20+90*(n_lam-1)))*1e6;

params.fbg.N_s = N_s;
params.fbg.D = 400:20:(400+20*(N_s-1));
params.fbg.periods = 530.73e-9 * ones(1, N_s);
params.fbg.deltaneffs = deltaneff * ones(1, N_s);
params.fbg.grating_lengths = 3e-3 * ones(1, N_s);

% Spectral spread: symmetric around center
params.fbg.lambda_shifts = linspace(-spread_nm/2, spread_nm/2, N_s) * 1e-9;

% Temperature gradient (additive)
if grad_nm > 0
    params.fbg.lambda_shifts = params.fbg.lambda_shifts + linspace(0, grad_nm, N_s) * 1e-9;
end

try
    out = runSimulation(params);
    mae = out.MAE;
    if isnan(mae); mae = 1e4; end
catch
    mae = 1e4;
end

end
