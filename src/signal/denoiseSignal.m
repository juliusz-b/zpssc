function [out] = denoiseSignal(x, denoise_m)
% DENOISESIGNAL Apply denoising to signal using specified method
%
% [out] = denoiseSignal(x, denoise_m)
%
% Inputs:
%   x - Input signal to denoise
%   denoise_m - Structure with denoising parameters:
%       denoise_m.method - 'wavelet', 'tvd', 'sgolay', or 'none'
%       Additional method-specific parameters:
%         - tvd: lambda, iterations
%         - sgolay: order, framelen
%
% Outputs:
%   out - Denoised signal

% Initialize output
out = [];

% Apply selected denoising method
switch denoise_m.method
    case 'wavelet'
        % Wavelet denoising
        for i = 1:size(x, 1)
            out(i, :) = wdenoise(x(i, :));
        end
        
    case 'tvd'
        % Total Variation Denoising
        for i = 1:size(x, 1)
            out(i, :) = tvd_mm(x(i, :), denoise_m.lambda, denoise_m.iterations);
        end
        
    case 'sgolay'
        % Savitzky-Golay filtering
        for i = 1:size(x, 1)
            out(i, :) = sgolayfilt(x(i, :), denoise_m.order, denoise_m.framelen);
        end
        
    case 'none'
        % No denoising, pass through
        out = x;
        
    otherwise
        error('Unsupported denoising method. Use: wavelet, tvd, sgolay, or none');
end
end
