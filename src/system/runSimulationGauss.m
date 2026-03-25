function [out] = runSimulationGauss(params)
% RUNSIMULATIONGAUSS Run FBG CDM simulation with Gaussian fit for λ_B estimation
%
% Identyczna z runSimulation() ale używa Gaussian fit + spline interpolation
% zamiast centroidu do estymacji λ_B.
%
% Różnica vs runSimulation:
%   - Centroid: λ_B = Σ R(λ)·λ / Σ R(λ) - wrażliwy na sidelobes (floor)
%   - Gaussian: fit A·exp(−(λ−μ)²/2σ²) - odporny na floor
%
% Użycie: identyczne jak runSimulation(params)
% Porównanie: uruchom oba i porównaj out.MAE

% Uruchom standardową symulację (cały pipeline)
out = runSimulation(params);

% Nadpisz estymację λ_B: Gaussian fit zamiast centroidu
neff = 1.447;
lambdas = out.lambdas;
lambdas_nm = lambdas * 1e9;
N_s = params.fbg.N_s;

lB_gauss = zeros(1, N_s);

for i = 1:N_s
    y = out.R_received(i, :);

    % Spline interpolation 10×
    x_interp = linspace(lambdas_nm(1), lambdas_nm(end), length(lambdas_nm) * 10);
    y_interp = interp1(lambdas_nm, y, x_interp, 'spline');
    y_interp = max(0, y_interp); % clamp negative

    % Threshold: 30% of max (jak w centroidzie)
    mx = max(y_interp) * 0.3;
    I = y_interp >= mx;
    x_fit = x_interp(I);
    y_fit = y_interp(I);

    if length(x_fit) >= 3
        % Gaussian model: A · exp(−(x−μ)² / (2σ²))
        gauss_fn = @(p, x) p(1) * exp(-(x - p(2)).^2 / (2 * p(3)^2));

        % Initial guess
        mu0 = sum(y_fit .* x_fit) / sum(y_fit); % centroid as starting point
        A0 = max(y_fit);
        sigma0 = 0.05; % ~50 pm

        try
            opts = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxIterations', 100);
            p_fit = lsqcurvefit(gauss_fn, [A0, mu0, sigma0], x_fit, y_fit, ...
                [0, lambdas_nm(1), 0.001], [Inf, lambdas_nm(end), 1], opts);
            lB_gauss(i) = p_fit(2);
        catch
            % Fallback: centroid
            lB_gauss(i) = mu0;
        end
    else
        % Too few points: use centroid
        mx2 = max(y) * 0.3;
        I2 = y >= mx2;
        if any(I2)
            lB_gauss(i) = sum(y(I2) .* lambdas_nm(I2)) / sum(y(I2));
        else
            lB_gauss(i) = NaN;
        end
    end
end

% Overwrite λ_B results with Gaussian fit
out.lB_rec_centroid = out.lB_rec;    % zachowaj centroid
out.MAE_centroid = out.MAE;          % zachowaj MAE centroid
out.lB_rec = lB_gauss;               % nadpisz Gaussian
out.lB_errors = lB_gauss - out.lB_org;
out.MAE = mean(abs(out.lB_errors), 'omitnan') * 1000;
out.MAX_err = max(abs(out.lB_errors)) * 1000;
out.estimation_method = 'gaussian_fit';

end
