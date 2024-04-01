function [filter_data] = generateInitialFilter(initial_filter, Fsample, Nb, L)
% GENERATEINITIALFILTER Generate initial filter for signal conditioning
%
% [filter_data] = generateInitialFilter(initial_filter, Fsample, Nb, L)
%
% Inputs:
%   initial_filter - Structure with filter parameters:
%       initial_filter.type - 'rrcos', 'matched', or 'none'
%       initial_filter.alpha - Roll-off factor for RRCOS filter
%   Fsample - Sampling frequency [Hz]
%   Nb - Samples per bit
%   L - Filter length in symbols
%
% Outputs:
%   filter_data - Filter coefficients

switch initial_filter.type
    case 'rrcos'
        % Root-raised cosine filter
        filter_data = rcosdesign(initial_filter.alpha, L, Nb);
        
    case 'matched'
        % Matched filter (simple rectangular pulse)
        filter_data = ones(1, Nb) / Nb;
        
    case 'none'
        % No filtering
        filter_data = 1;
        
    otherwise
        error('Unsupported filter type. Use: rrcos, matched, or none');
end
end
