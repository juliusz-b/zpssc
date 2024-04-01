function [out] = myPrbs(u, poly, initial_poly)
% MYPRBS Generate pseudo-random binary sequence
%
% [out] = myPrbs(u, poly, initial_poly)
%
% Inputs:
%   u - Order of the sequence (resulting in sequence of length 2^u-1)
%   poly - (optional) Polynomial coefficients for LFSR
%   initial_poly - (optional) Initial state of LFSR
%
% Outputs:
%   out - Generated PRBS sequence

% Check valid u values
if u < 4 || (u > 23 && (u ~= 21 && u ~= 22))
    error('PRBS order must be between 4 and 23 (except 21 and 22)');
end

% Set default polynomial and initial state if not provided
if nargin < 2
    switch u
        case 4
            poly = [1 1 0 0 1];
        case 5
            poly = [1 0 1 0 0 1];
        case 6
            poly = [1 1 0 0 0 0 1];
        case 7
            poly = [1 1 0 0 0 0 0 1];
        case 8
            poly = oct2poly(435);  % [1 0 1 1 1 0 0 0 1]
        case 9
            poly = [1 0 0 0 1 0 0 0 0 1];
        case 10
            poly = [1 0 0 1 0 0 0 0 0 0 1];
        case 11
            poly = [1 0 1 0 0 0 0 0 0 0 0 1];
        case 12
            poly = [1 1 1 0 0 0 0 0 1 0 0 0 1];
        case 13
            poly = [1 1 1 0 0 1 0 0 0 0 0 0 0 1];
        case 14
            poly = [1 1 1 0 0 0 0 0 0 0 0 0 1 0 1];
        case 15
            poly = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
        case 16
            poly = [1 0 1 1 0 1 0 0 0 0 0 0 0 0 0 0 1];
        case 17
            poly = [1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
        case 18
            poly = [1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1];
        case 19
            poly = [1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
        case 20
            poly = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1];
        case 23
            poly = [1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
    end
    
    % Default initial state - all zeros with a single 1 at the end
    initial_poly = [zeros(1, u-1) 1];
else
    % If poly provided but no initial state
    if nargin < 3
        if length(poly) ~= u+1
            error('Wrong length of polynomial for given order!');
        end
        % Reverse the polynomial coefficients for initial state
        initial_poly = fliplr(poly(1:end-1));
    else
        % Check if initial state has correct length
        if length(poly)-1 ~= length(initial_poly)
            error('Wrong length of initial_poly!');
        end
    end
end

% Generate PRBS sequence using LFSR
out = lfsr(poly, initial_poly, 2^u-1);
end

function [out_code] = lfsr(polynomial, beginning_polynomial, number_of_iterations)
% LFSR - Linear Feedback Shift Register implementation
%
% Implements LFSR with given polynomial coefficients and initial state

if ~isvector(polynomial)
    error('polynomial must be a vector');
end

L = length(polynomial)-1;
if length(beginning_polynomial) ~= L
    error('beginning_polynomial must have same length as polynomial-1');
end

reg = beginning_polynomial;

flipped_polynomial = flip(polynomial(1:end-1));
flipped_polynomial_temp = flipped_polynomial;
ones_L = sum(flipped_polynomial);

out_code = zeros(1, number_of_iterations);
for it = 1:number_of_iterations
    flipped_polynomial_temp = flipped_polynomial;
    
    if ones_L > 1
        % Multiple feedback taps
        for j = 1:ones_L-1
            if j == 1
                % Find first tap
                id1 = find(flip(flipped_polynomial_temp) == 1, 1, 'first');
                flipped_polynomial_temp = flip(flipped_polynomial_temp);
                flipped_polynomial_temp(id1) = 0;
                flipped_polynomial_temp = flip(flipped_polynomial_temp);
                id1 = L-id1+1;
                
                % Find second tap
                id2 = find(flip(flipped_polynomial_temp) == 1, 1, 'first');
                flipped_polynomial_temp = flip(flipped_polynomial_temp);
                flipped_polynomial_temp(id2) = 0;
                flipped_polynomial_temp = flip(flipped_polynomial_temp);
                id2 = L-id2+1;
                
                % XOR the first two taps
                bit = xor(reg(id1), reg(id2));
            else
                % Find next tap and XOR with current bit
                id1 = find(flip(flipped_polynomial_temp) == 1, 1, 'first');
                flipped_polynomial_temp = flip(flipped_polynomial_temp);
                flipped_polynomial_temp(id1) = 0;
                flipped_polynomial_temp = flip(flipped_polynomial_temp);
                id1 = L-id1+1;
                bit = xor(bit, reg(id1));
            end
        end
        
        % Output current register value
        out_code(it) = reg(end);
        
        % Shift register and insert new bit
        reg = circshift(reg, 1);
        reg(1) = bit;
    else
        % Single feedback tap
        out_code(it) = reg(end);
        reg = circshift(reg, 1);
        
        if sum(flipped_polynomial == 1) == 0
            reg(1) = reg(end);
        else
            reg(1) = reg(flipped_polynomial == 1);
        end
    end
end
end

function poly = oct2poly(oct_val)
% Convert octal number to polynomial representation
bin_str = dec2bin(base2dec(num2str(oct_val), 8));
poly = zeros(1, length(bin_str) + 1);
for i = 1:length(bin_str)
    poly(i) = str2double(bin_str(i));
end
poly(end) = 1;  % Add the x^0 term
end
