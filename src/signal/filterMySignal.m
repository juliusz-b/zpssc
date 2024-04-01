function [data_out] = filterMySignal(data_in, filter_data)
% FILTERMYSIGNAL Apply filtering to signal data
%
% [data_out] = filterMySignal(data_in, filter_data)
%
% Inputs:
%   data_in - Input signal data
%   filter_data - Filter coefficients
%
% Outputs:
%   data_out - Filtered signal data

% Initialize output
data_out = [];

% Process each row of input data
for i = 1:size(data_in, 1)
    % Apply filter via convolution
    flt_data = conv(data_in(i, :), filter_data);
    
    % Find delay introduced by filter
    D = finddelay(data_in(i, :), flt_data) + 1;
    
    % Find delay from the right side
    D2 = finddelay(fliplr(data_in(i, :)), fliplr(flt_data(1:end)));
    
    % Extract filtered data with delay compensation
    flt_data = flt_data(D:end-D2);
    
    % Ensure output length matches input
    if length(flt_data) > size(data_in, 2)
        flt_data = flt_data(1:end-length(flt_data)+size(data_in, 2));
    end
    
    % Add to output
    data_out = [data_out; flt_data]; 
end
end
