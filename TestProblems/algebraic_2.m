% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

function y_eval = algebraic_2(theta)
% ALGEBRAIC_2: Function evaluates f(x) = (1 + sum(c_k*xi_k))^(-(d+1)) 
% for a set of input theta parameters in [-1,1]^d

% Determine the number of samples and parameters 
[N,d] = size(theta);

% Declare avector for the function evaluations
y_eval = zeros(N,1);

% Evaluate the function for all theta parameters
for i=1:N
    
    % Evaluate the component coefficients
    c_k_comp = 1./((1:d).^2);
    c_k_vect = c_k_comp/(4*sum(c_k_comp));
        
    % Transform inputs to [0,1] and compute inner product with c_k
    new_theta = 0.5*theta(i,:) + 0.5;
    algebraic_fn = sum(c_k_vect.*new_theta);
    
    % Assign result to y_eval
    y_eval(i) = (1 + algebraic_fn)^(-1*(d+1));

end

% -- END OF FILE --
