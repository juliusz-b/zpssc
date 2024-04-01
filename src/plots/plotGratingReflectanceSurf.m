function plotGratingReflectanceSurf(lambdas, R_original, R_received, length_domain, dane_po_przetworzeniu, dane_w_kanale, dane_referencyjne, N_s, D_s, R_simulated, simulation_range, new_lmb, R_filtered)
% PLOTGRATINGREFLECTANCESURF Visualize FBG reflectance as a surface plot
%
% plotGratingReflectanceSurf(lambdas, R_original, R_received, length_domain, 
%                          dane_po_przetworzeniu, dane_w_kanale, dane_referencyjne,
%                          N_s, D_s, R_simulated, simulation_range, new_lmb, R_filtered)
%
% Inputs:
%   lambdas - Wavelength array [m]
%   R_original - Original FBG reflectance spectra
%   R_received - Detected FBG reflectance spectra
%   length_domain - Distance domain [m]
%   dane_po_przetworzeniu - Processed data
%   dane_w_kanale - Channel data
%   dane_referencyjne - Reference data
%   N_s - Number of gratings
%   D_s - Grating positions
%   R_simulated - Simulated high-resolution reflectance spectra
%   simulation_range - Wavelength range for simulation [m]
%   new_lmb - Upsampled wavelength array [m]
%   R_filtered - Filtered reflectance spectra

% Create 3D visualization of data
Z = zeros(size(dane_po_przetworzeniu, 1), length(length_domain));
for i = 1:size(dane_po_przetworzeniu, 1)
    % Process data with variance filter
    xcv = varFilter(dane_po_przetworzeniu(i, :), size(dane_referencyjne, 2), mean(diff(D_s)), true, true, false);
    
    % Apply filter to focus on grating positions
    xcv_filter = ones(size(xcv));
    xcv_filter(1:D_s(1)-1) = 0;
    xcv = xcv .* xcv_filter;
    
    % Store in output matrix
    Z(i, :) = xcv;
end

% Create a matrix without envelope effects
Z3 = dane_po_przetworzeniu - movmean(dane_po_przetworzeniu, D_s(2) - D_s(1), 2);

% Create figure for 3D plot
figure('Renderer', 'opengl', 'units', 'normalized', 'Color', 'w');

% Plot markers at grating positions
for i = 1:length(length_domain)
    if length_domain(i) >= 390 && length_domain(i) < 530
        p3 = plot3(repmat(length_domain(i), length(lambdas), 1), lambdas*1e9, Z3(:, i), '*', 'linewidth', 2, 'color', 'black');
        zlim([0.2 2]);
        hold on;
    end
end

% Create surface plot of correlation data
figure('Renderer', 'opengl', 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', 'w');
surf(lambdas*1e9, length_domain, Z', 'EdgeColor', 'none');
ylim([390 530]);
view([0 90]);
shading interp;
colormap jet;
xlabel('Wavelength [nm]');
ylabel('Position [m]');
title('Cross-correlation Surface');
colorbar;