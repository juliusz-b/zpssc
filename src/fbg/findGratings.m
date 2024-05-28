function [delay, R] = findGratings(xcv, reference_data, length_domain, N_s, D_s)
% FINDGRATINGS Locate FBG gratings from correlation data
%
% [delay, R] = findGratings(xcv, reference_data, length_domain, N_s, D_s)
%
% Inputs:
%   xcv - Cross-correlation data
%   reference_data - Original reference signals
%   length_domain - Distance domain
%   N_s - Number of gratings
%   D_s - Expected grating positions
%
% Outputs:
%   delay - Detected delay for each grating
%   R - Reflectance values for each grating and wavelength

% Number of gratings to detect
N_of_gratings = N_s;

% Initialize output arrays
R = zeros(N_of_gratings, size(xcv, 1));
delay = zeros(1, N_of_gratings);

% Set detection mode
MODE = 1;
max_peaks = -Inf;

% Process each wavelength
for i = 1:size(xcv, 1)
    % Process correlation data for current wavelength
    %xcv_filtered = varFilter(xcv(i,:), size(reference_data, 2), mean(diff(D_s))*1.5, true, true, true);
    xcv_filtered = xcv(i,:);
    
    % Non-zero indices that are within the expected grating locations
    NZIX = xcv_filtered ~= 0;
    xcv_filtered(NZIX) = xcv_filtered(NZIX) ./ xcv_filtered(NZIX) .* xcv(i, NZIX);
    
    % Find peaks in the filtered data
    [peaks_locs, peaks_, ~] = findPeak(xcv_filtered, N_of_gratings);
    
    % Calculate reference peak value
    if sum(reference_data) == length(reference_data)
        peaks_locs_ref = 1;
        peaks_ref = 1;
    else
        [peaks_locs_ref, peaks_ref, ~] = findPeak(calcXcov(reference_data(i,:)), 1,1);
    end
    
    % Normalize peaks by reference peak
    peaks_ = peaks_ ./ peaks_ref;
    
    % Store reflectances for this wavelength
    R(:, i) = peaks_;
    
    % Determine delay based on detection mode
    if MODE == 1
        if sum(peaks_) > max_peaks
            max_peaks = sum(peaks_);
            delay = length_domain(peaks_locs);
        end
    else
        if i == 1
            delay = length_domain(peaks_locs);
        else
            delay = (delay + length_domain(peaks_locs)) / 2;
        end
    end
end

% Ensure reflectance values are positive
R = abs(R);
end

function [peaks_locs, peaks_, sort_idx] = findPeak(x, N_of_gratings, D_s)
% Find peaks in signal x

if nargin < 3
    MODE = 'normal';
    if nargin < 2
       N_of_gratings = 1; 
    end
else
   MODE = 'reference';
   if length(D_s) < 2
    sample_diff = (D_s) / 2;
   else
    sample_diff = round(mean(diff(D_s)) / 2);
   end
end

switch MODE
    case 'normal'
        % Shift by 10 samples to ensure we catch the first sample
        SHIFT = 10;
        xx = 1:length(x);
        xx_cs = circshift(xx, SHIFT);
        
        % Find peaks
        [peaks_, peaks_locs] = findpeaks(circshift(x, SHIFT), xx);
        
        % Remove peaks before the first expected location
        %peaks_(peaks_locs < D_s(1)) = [];
        %peaks_locs(peaks_locs < D_s(1)) = [];
        
        % Sort by amplitude
        [peaks_, sort_idx] = sort(peaks_, 'descend');
        peaks_locs = peaks_locs(sort_idx);

        % Take the top N peaks
        peaks_ = peaks_(1:min(N_of_gratings, length(peaks_)));
        peaks_locs = peaks_locs(1:min(N_of_gratings, length(peaks_locs)));

        % Sort by location
        [peaks_locs, sort_idx] = sort(peaks_locs);
        peaks_ = peaks_(sort_idx);

        % Convert to original indexing
        peaks_locs = xx_cs(peaks_locs);
        
    case 'reference'
        % For reference mode, search near expected positions
        xx = 1:length(x);
        peaks_ = zeros(1, length(D_s));
        peaks_locs = zeros(1, length(D_s));
        
        for i = 1:length(D_s)
            d = D_s(i);
            
            % Define search window
            if i == 1
                ind = 1:d+sample_diff;
            else
                ind = D_s(i-1)+sample_diff:d+sample_diff;
            end
            
            % Ensure indices are within valid range
            xx_it = xx(ind(ind <= length(xx)));
            
            % Find maximum in window
            [mx, it] = max(x(xx_it));
            peaks_(i) = mx;
            peaks_locs(i) = xx_it(it);
        end
        
        sort_idx = 0;
end
end
