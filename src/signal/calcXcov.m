function [xcov_out] = calcXcov(data_in, data_original, should_ref, xcov_block)
% CALCXCOV Calculate cross-covariance between signals
%
% [xcov_out] = calcXcov(data_in, data_original, should_ref, xcov_block)
%
% Inputs:
%   data_in - Input data signal
%   data_original - Reference signal
%   should_ref - (optional) Boolean flag for reference mode
%   xcov_block - (optional) Parameters for block-based processing
%
% Outputs:
%   xcov_out - Cross-covariance result

% Default parameters
if nargin < 4
    xcov_block = {'normal', 1, 1};
    if nargin < 3
        should_ref = false;
        if nargin < 2
           data_original = data_in; 
        end
    end
end

% Check if we can use the standard mode
if size(data_in, 1) == size(data_original, 1)
   should_ref = false; 
end

% Get processing mode
xcov_mode = xcov_block{1};
if ~contains(xcov_mode, {'mean', 'normal'})
    error('Unsupported mode: use "mean" or "normal"');
end

% Determine the size of the output
M = max(size(data_in, 2), size(data_original, 2));
xcov_out = zeros(size(data_original, 1), M);

% Check input dimensions
if (size(data_in) ~= size(data_original))
    if ~should_ref
        error('Data dimensions mismatch when not in reference mode');
    end
end

% Number of reference signals
s_ = size(data_original, 1);

% Ensure data_original has sufficient length
if size(data_original, 2) < size(data_in, 2)
    data_original_padded = [data_original, zeros(size(data_original, 1), size(data_in, 2) - size(data_original, 2))];
else
    data_original_padded = data_original;
end

% Calculate cross-covariance
for i = 1:s_
    if should_ref
        % Cross-covariance between input and reference signal
        xcov_temp = xcorr(data_in(1, :), data_original_padded(i, :));
        xcov_out(i, :) = xcov_temp(floor(length(xcov_temp)/2)+1:end);
    else
        % Cross-covariance between corresponding signals
        xcov_temp = xcorr(data_in(i, :), data_original_padded(i, :));
        xcov_out(i, :) = xcov_temp(floor(length(xcov_temp)/2)+1:end);
    end
end
end
