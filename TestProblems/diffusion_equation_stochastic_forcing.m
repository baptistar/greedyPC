% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

function soln_x_loc = diffusion_equation_stochastic_forcing(theta)
% DIFFUSION_EQUATION_STOCHASTIC_FORCING: Function evaluates the 2D 
% stochastic diffusion equation with KL expansions for the diffusivity
% and the forcing term at different theta values based on the created 
% model and the declared parameters.

% flip theta
theta = theta';

%% Declare Constants
 
% Constants for the Diffusion Equation
sigma = 0.7;
corl  = [0.05, 0.05];

% Quantity of Interest in Spatial Domain
x_loc = [0.5, 0.5];

% Define 1D Spatial Dimensions
dim_lim = [0,1];

% Determine the number of parameters and points
n_par = size(theta,1);
n_pts = size(theta,2);

%% Evaluate Governing Equations

% Find file_name
file_name = ['../TestProblems/Diffusion_SForce_s' num2str(sigma) '_p' num2str(n_par) '.mat'];

% Extract global K and F matrices and mesh terms (p & t)
if exist(file_name , 'file')
    load(file_name)
else
    [~,global_K,global_F,~,~,~,~,int_p] = diff_2D_stochastic_forcing(n_par/2, sigma, corl, dim_lim);
    save(file_name,'global_K','global_F','int_p');
end

%% Evaluate Diffusion PDE

% Declare vector for solution
soln_x_loc = zeros(n_pts,1);

% Divide theta parameters
theta_M = theta(1:n_par/2,:);
theta_F = theta(n_par/2+1:end,:);

% Extract triangular cell containing interpolation point
X = int_p(1,:)';
Y = int_p(2,:)';
tri = delaunay (X, Y);
idx = tsearchn([X, Y], tri, [x_loc(1), x_loc(2)]);
pts = tri(idx, :);

for i=1:n_pts
    
    % Evaluate the leading matrix
    K_matrix = global_K{1};
    for j=1:n_par/2
        K_matrix = K_matrix + global_K{j+1}*theta_M(j,i);
    end
    
    % Evaluate the source term
    F_matrix = global_F{1};
    for j=1:n_par/2
        F_matrix = F_matrix + global_F{j+1}*theta_F(j,i);
    end
    
    % Determine the exact solution
    soln = K_matrix\F_matrix;
    
    % Extract value at location of interest (x_loc) by interpolation
    F = scatteredInterpolant(X(pts), Y(pts), soln(pts));
    soln_x_loc(i) = F(x_loc(1), x_loc(2));
    
end

end

% -- END OF FILE --
