function [out] = varFilter(x, var_K, center_K, should_use_parfor, should_detrend, should_center)
% VARFILTER Variance-based filtering for signal enhancement
%
% [out] = varFilter(x, var_K, center_K, should_use_parfor, should_detrend, should_center)
%
% Inputs:
%   x - Input signal to filter
%   var_K - Window size for variance calculation
%   center_K - (optional) Window size for centering operation
%   should_use_parfor - (optional) Boolean flag to use parallel processing
%   should_detrend - (optional) Boolean flag to detrend signal
%   should_center - (optional) Boolean flag to center signal
%
% Outputs:
%   out - Filtered signal

% Initialize output
out = x;

% Set default parameters
if nargin < 6
    should_center = true;
    if nargin < 5
        should_detrend = true;
        if nargin < 4
            should_use_parfor = true;
            if nargin < 3
                center_K = var_K;
            end
        end
    end
end

% Configure parallel processing
if should_use_parfor
    workers = Inf;
else
    workers = 0;
end

% Process signal point by point
parfor (i = 1:length(x), workers)
    % Define window bounds
    L = max(1, i - round(var_K/2));
    R = min(length(x), i + round(var_K/2) - 1);
    
    % Extract data window
    data = x(L:R);
    
    % Apply detrending if requested
    if should_detrend
        [data, median_] = detrendStep(data, x(i), R, L);
    else
        median_ = x(i);
    end
    
    % Calculate variance-based metric
    out(i) = varStep(data, median_);
    
    % Apply centering if requested
    if should_center
        % Define centering window bounds
        L_cent = max(1, i - round(center_K/2));
        R_cent = min(length(x), i + round(center_K/2) - 1);
        
        % Extract centering window data
        if should_detrend
            [data_cent] = detrendStep(x(L_cent:R_cent), x(i), R_cent, L_cent);
        else
            data_cent = x(L_cent:R_cent);
        end
        
        % Find maximum in centered window
        [mx, it_mx] = max(data_cent);
        if it_mx == floor((R_cent-L_cent+1)/2)
            out(i) = mx;
        else
            out(i) = 0;
        end
    end
end
end

function out = varStep(data, median_)
% Calculate variance-based metric
denom = mean(data.^2);
nom = median_^2;
out = nom/denom;
end

function [data_out, median_out] = detrendStep(data, center_val, R, L)
% Detrend data window

% Calculate linear trend
alpha = mean(diff(data));
beta = mean(data);

% Remove linear trend
data_out = data - (0:length(data)-1) * alpha - beta;

% Apply additional detrending
data_out = detrend(data_out, 2);

% Calculate detrended center value
median_out = center_val - beta - (R-L)/2 * alpha;

% Use center of detrended data
median_out = data_out(ceil(end/2));
end
