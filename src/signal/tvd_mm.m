function [x, cost] = tvd_mm(y, lam, Nit)
% TVD_MM Total variation denoising using majorization-minimization
%
% [x, cost] = tvd_mm(y, lam, Nit)
%
% Total variation denoising of signal y using majorization-minimization
% algorithm and banded linear systems.
%
% Inputs:
%   y - Noisy signal
%   lam - Regularization parameter
%   Nit - Number of iterations
%
% Outputs:
%   x - Denoised signal
%   cost - Cost function history (optional)
%
% Reference:
% 'On total-variation denoising: A new majorization-minimization
% algorithm and an experimental comparison with wavalet denoising.'
% M. Figueiredo, J. Bioucas-Dias, J. P. Oliveira, and R. D. Nowak.
% Proc. IEEE Int. Conf. Image Processing, 2006.

% Ensure y is a column vector
y = y(:);

% Initialize cost function history
cost = zeros(1, Nit);

% Get signal length
N = length(y);

% Create sparse identity and difference matrices
I = speye(N);
D = I(2:N, :) - I(1:N-1, :);
DDT = D * D';

% Initialize solution with noisy signal
x = y;
Dx = D*x;
Dy = D*y;

% Iteratively solve the problem
for k = 1:Nit
    % Compute weight matrix based on current solution
    F = sparse(1:N-1, 1:N-1, abs(Dx)/lam) + DDT;
    
    % Solve the banded linear system
    x = y - D'*(F\Dy);
    
    % Update Dx for next iteration
    Dx = D*x;
    
    % Compute cost function value
    cost(k) = 0.5*sum(abs(x-y).^2) + lam*sum(abs(Dx));
end

% Return x as a row vector
x = x(:)';
end
