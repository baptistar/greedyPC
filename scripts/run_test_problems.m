% Polynomial Chaos Randomized Greedy Algorithm
%
% Copyright (C) 2019 R. Baptista & P. Nair
%
% GreedyPC is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% GreedyPC is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License 
% along with GreedyPC. If not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) 2019 MIT & University of Toronto
% Authors: Ricardo Baptista & Prasanth Nair
% E-mails: rsb@mit.edu or ricarsb@gmail.com & pbn@utias.utoronto.ca
%

% Reproduce studies in JCP Paper
clear; close all; clc

addpath('../SpectralToolbox')
addpath('../TestProblems')
addpath('../Methods')
addpath('../tools')

% make folder to store results
if ~exist('../results', 'dir')
       mkdir('../results')
end

%% ALGEBRAIC 1

% define test parameters
test = struct;
test.d          = 10;
test.order      = 4;
test.grid_level = 5;
test.basis      = 'Legendre';
test.func       = @(x) algebraic_1(x);
test.func_str   = 'algebraic_1';
test.N          = [50, 100, 200, 400, 600, 800, 1000];
test.n_trials   = 5;
test.MCtrials   = 100;

% run algorithms
run_sampling(test);
run_algorithms(test)
create_plots(test)

%% ALGEBRAIC 2

% define test parameters
test = struct;
test.d          = 10;
test.order      = 4;
test.grid_level = 5;
test.basis      = 'Legendre';
test.func       = @(x) algebraic_2(x);
test.func_str   = 'algebraic_2';
test.N          = [50, 100, 200, 400, 600, 800, 1000];
test.n_trials   = 5;
test.MCtrials   = 100;

% run algorithms
run_sampling(test);
run_algorithms(test)
create_plots(test)

%% ALGEBRAIC 3

% define test parameters
test = struct;
test.d          = 10;
test.order      = 4;
test.grid_level = 5;
test.basis      = 'Legendre';
test.func       = @(x) algebraic_3(x);
test.func_str   = 'algebraic_3';
test.N          = [50, 100, 200, 400, 600, 800, 1000];
test.n_trials   = 5;
test.MCtrials   = 100;

% run algorithms
run_sampling(test);
run_algorithms(test)
create_plots(test)

%% DIFFUSION EQUATION

% define test parameters
test = struct;
test.d          = 20;
test.order      = 3;
test.grid_level = 5;
test.basis      = 'Legendre';
test.func       = @(x) diffusion_equation_stochastic_forcing(x);
test.func_str   = 'diff_equation';
test.N          = [50, 100, 200, 400, 600, 800, 1000];
test.n_trials   = 5;
test.MCtrials   = 100;

% run algorithms
run_sampling(test);
run_algorithms_lowmem(test)
create_plots(test)

% -- END OF FILE --
