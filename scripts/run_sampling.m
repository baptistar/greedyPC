% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

function run_sampling(test)
% RUN_SAMPLING: Runs the sampling based methods (MC and QMC)
% on the test case specified in the test struct. The results 
% are saved in the results/test_name file.

% extract parameters
d          = test.d;
order      = test.order;
grid_level = test.grid_level;
basis      = test.basis;
func       = test.func;
func_str   = test.func_str;
N          = test.N;
n_trials   = test.MCtrials;

% call sparse grid to determine the test samples
[XTest, wTest] = sparse_grid(d, grid_level, 'unif');
YTest = func(XTest);

% evaluate mean and standard deviation of yTest
meanTest = wTest*YTest;
stdTest  = sqrt(wTest*YTest.^2 - meanTest^2);

% declare vectors to store results
MC_out  = cell(length(N), n_trials);
QMC_out = cell(length(N), 1);

% define grid of parameters
[Run, Nm] = meshgrid(1:n_trials, N);

% Run for each algorithm
n_proc = 10;
c = parcluster('local');
c.NumWorkers = n_proc;
parpool(n_proc);

%% Monte Carlo Sampling
parfor i=1:length(N)*n_trials
    
    % setup structs to save results
    MCpp  = struct;
    
    % generate training data: uniform inputs on [-1,1]
    X = 2*rand(Nm(i), d) - 1;
    Y = func(X);

    % evaluate statistics
    meanMC = mean(Y);
    stdMC  = sqrt(mean(Y.^2) - meanMC^2);

    % Evaluate statistics
    MCpp.MeanE = norm(meanMC - meanTest)/norm(meanTest);
    MCpp.StdE  = norm(stdMC - stdTest)/norm(stdTest);
    
    % save struct
    MC_out{i} = MCpp;

end

% check if sobol is in current path
if exist('i4_sobol_generate','file') ~= 2
    error('Please add Sobol scripts to path - see README')
end

%% Quasi Monte Carlo (Sobol) Sampling
parfor i=1:length(N)

    % setup structs to save results
    QMCpp  = struct;
    
    % Generate samples
    X = 2*i4_sobol_generate(d, N(i), 0)' - 1;
    Y = func(X);
    
    % evaluate statistics
    meanQMC = mean(Y);
    stdQMC  = sqrt(mean(Y.^2) - meanQMC^2);

    % Evaluate statistics
    QMCpp.MeanE = norm(meanQMC - meanTest)/norm(meanTest);
    QMCpp.StdE  = norm(stdQMC - stdTest)/norm(stdTest);

    % save struct
    QMC_out{i} = QMCpp;

end

delete(gcp('nocreate'))

% save matlab file with results
save(['../results/' func_str '_d' num2str(d) '_samples']);

end
