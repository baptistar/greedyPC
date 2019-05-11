% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

function [nodes, weights] = sparse_grid(n_dim, grid_level, dist)
% SPARSE_GRI: Function computes the nodes and weights to evaluate a 
% multi-dimensional integral using a Cleanshaw-Curtis sparse grid scheme.
%
% This function uses the sparse_grid_cfn_size.m and sparse_grid_cc.m 
% functions that are supplied with the SPARSE_GRID_CC program
% by John Burkardt (2009) available at:
% people.sc.fsu.edu/~jburkardt/m_src/sparse_grid_cc/sparse_grid_cc.html

% Declare File Name
file_name = ['sparsegrid_cc_dim_' num2str(n_dim)  '_level_' num2str(grid_level) '.mat'];

% Check if file exists, if not run sparse grid code
if exist(file_name, 'file')

	% Load file contents
	load(file_name);
    
else

	% check if sparse grid tool is in path
	if (exist('sparse_grid_cfn_size','file') ~= 2) || ...
	   (exist('sparse_grid_cc','file') ~= 2)
		error('Please add sparse_grid tools to current path - see README');
	end

	% Calculate the number of points in the sparse grid
	num_pts = sparse_grid_cfn_size(n_dim, grid_level);

	% Generate the nodes and weights of the sparse grid quadrature rule
	[weights, nodes] = sparse_grid_cc(n_dim, grid_level, num_pts);

	% Save files
	save(file_name, 'weights', 'nodes');

end

nodes = nodes';

% Scale weights based on distribution
if strcmp(dist,'unif')
	weights = weights/(2^n_dim);
else
	error('Distribution is Not Recognized')
end