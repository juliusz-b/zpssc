function [out] = myGold(m)
% MYGOLD Generate Gold sequence set
%
% [out] = myGold(m)
%
% Inputs:
%   m - Order of the sequence (resulting in sequences of length 2^m-1)
%
% Outputs:
%   out - Matrix containing the generated Gold sequences

% Check if m is a valid value for Gold sequences
valid_m = [4, 5, 6, 7, 8, 9, 10, 11, 15];
if isempty(find(valid_m == m, 1, 'first'))
    error('Invalid m value for Gold sequences. Valid values: 4,5,6,7,8,9,10,11,15');
end

% Select appropriate primitive polynomials based on m
switch m
    case 4
        poly1 = oct2poly(23);  % Octal representation of polynomial
    case 5
        poly1 = oct2poly(45);
    case 6
        poly1 = oct2poly(103);
    case 7
        poly1 = oct2poly(211);
    case 8
        poly1 = oct2poly(435);
    case 9
        poly1 = oct2poly(1021);
    case 10
        poly1 = oct2poly(2011);
    case 11
        poly1 = oct2poly(4005);
    case 15
        poly1 = oct2poly(100003);
end

% Generate Gold sequence using decimation method
temp = gold_generator(m, poly1, fliplr(poly1(1:end-1)));

% Sequence length
L = 2^m - 1;

% Extract individual sequences from the generated set
out = zeros(2^m + 1, L);
for i = 1:(2^m + 1)
    out(i, :) = temp((L)*(i-1)+1:(L)*i);
end
end

function [out] = gold_generator(m, poly1, initial_state)
% Generate Gold sequences using preferred pair of m-sequences

% Sequence length
L = 2^m - 1;

% Generate first m-sequence
out1 = lfsr_generator(poly1, initial_state, L);

% For Gold sequences, we need a second m-sequence with appropriate decimation
% The decimation value should ensure the two sequences have a three-valued
% cross-correlation function
if m == 4
    k = 3;
elseif m == 8
    k = 4;
else
    % Find appropriate decimation value
    for k = 1:20
        if gcd(m, k) == 1 && mod(m, 2) == 1
            break;
        elseif gcd(m, k) == 2 && mod(m, 4) == 2
            break;
        end
    end
end

% Generate decimated sequence
q = 2^k + 1;
out2 = zeros(1, L);
for i = 1:L
    % Decimation by q
    idx = mod((i-1)*q, L) + 1;
    out2(i) = out1(idx);
end

% Generate all Gold sequences by XORing the first sequence with all cyclic shifts of the second
out = zeros(1, (2^m + 1) * L);

% First two sequences are the original m-sequences
out(1:L) = out1;
out(L+1:2*L) = out2;

% Generate Gold sequences by XORing with cyclic shifts
for i = 1:(2^m - 1)
    shifted = circshift(out2, i-1);
    out(2*L+1+(i-1)*L:2*L+i*L) = xor(out1, shifted);
end
end

function [out] = lfsr_generator(polynomial, initial_state, length_out)
% Generate sequence from LFSR with given polynomial and initial state

% Length of the shift register
reg_length = length(polynomial) - 1;

% Initialize register with initial state
register = initial_state;

% Generate sequence
out = zeros(1, length_out);
for i = 1:length_out
    % Output last bit of register
    out(i) = register(end);
    
    % Calculate feedback bit
    feedback = 0;
    for j = 1:reg_length
        if polynomial(j) == 1
            feedback = xor(feedback, register(j));
        end
    end
    
    % Shift register and add feedback
    register = [feedback, register(1:end-1)];
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
