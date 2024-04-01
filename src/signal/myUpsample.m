function [data_out] = myUpsample(data, upsample_ratio)
% MYUPSAMPLE Upsample data by repeating each sample
%
% [data_out] = myUpsample(data, upsample_ratio)
%
% Inputs:
%   data - Input data matrix or vector
%   upsample_ratio - Number of times to repeat each sample
%
% Outputs:
%   data_out - Upsampled data

% Initialize output data
data_out = zeros(size(data, 1), size(data, 2) * upsample_ratio);

% Upsample each row of data
for i = 1:size(data, 1)
    for j = 1:size(data, 2)
        % Repeat each sample upsample_ratio times
        data_out(i, upsample_ratio * (j - 1) + 1:j * upsample_ratio) = data(i, j);
    end
end
end
