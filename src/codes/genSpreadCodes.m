function [codes_out] = genSpreadCodes(L, U, code_type)
% GENSPREADCODES Generate spreading codes of various types
%
% [codes_out] = genSpreadCodes(L, U, code_type)
%
% Inputs:
%   L - Length of the code
%   U - Number of users/codes to generate
%   code_type - Type of code to generate:
%       'gold' - Gold sequences
%       'walsh' - Walsh-Hadamard codes
%       'kasami' - Kasami sequences
%       'prbs' - Pseudo-Random Binary Sequences
%       'randi' - Random binary sequences
%
% Outputs:
%   codes_out - Matrix of generated codes (U x L)

% Set random seed for reproducibility
rng('default');

% Initialize output
codes_out = [];

% Convert code type to lowercase
code_type = lower(code_type);

% Generate the appropriate type of code
switch code_type
    case 'gold'
        % Gold sequences using the myGold function
        seq = myGold(nextpow2(L));
        
    case 'walsh'
        % Walsh-Hadamard codes
        seq = hadamard(2^nextpow2(L*U));
        
    case 'kasami'
        % Kasami sequences using the myKasami function
        seq = myKasami(nextpow2(L));
        
    case 'prbs'
        % Pseudo-Random Binary Sequences
        tmp = myPrbs(nextpow2(L*U));
        
        % Distribute the PRBS sequence among users
        seq = zeros(U, L);
        for i = 1:U
            tmp2 = tmp(i:U:end);
            seq(i,:) = tmp2(1:L);
        end
        
        % Return directly since we've already formatted the output
        codes_out = seq;
        return;
        
    case 'randi'
        % Random binary sequences
        seq = randi([0 1], [L U])';
        
    otherwise
        error('Unknown code type. Supported types: gold, walsh, kasami, prbs, randi');
end

% Extract codes for each user
for i = 1:U
    if strcmp(code_type, 'gold')
        % For Gold sequences, skip the first entries which are the base sequences
        if i+6 <= size(seq, 1)
            codes_out(i,:) = seq(i+6, 1:L);
        else
            % If not enough codes, wrap around
            codes_out(i,:) = seq(mod(i+6-1, size(seq, 1))+1, 1:L);
        end
    else
        % For other code types, use sequential codes
        if i <= size(seq, 1)
            codes_out(i,:) = seq(i, 1:L);
        else
            % If not enough codes, wrap around
            codes_out(i,:) = seq(mod(i-1, size(seq, 1))+1, 1:L);
        end
    end
end
end
