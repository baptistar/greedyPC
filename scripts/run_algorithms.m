% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

function run_algorithms(test)
% RUN_ALGORITHMS: Runs the linear regression algorithms
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
n_trials   = test.n_trials;

% call sparse grid to determine the test samples
[XTest, wTest] = sparse_grid(d, grid_level, 'unif');
YTest = func(XTest);

% declare vectors to store results
FPC_out  = cell(length(N), n_trials);
BPDN_out = cell(length(N), n_trials);
OMP_out  = cell(length(N), n_trials);
OMPN_out = cell(length(N), n_trials);
RGA_out  = cell(length(N), n_trials);

% define grid of parameters
[Run, Nm] = meshgrid(1:n_trials, N);

% Run for each algorithm
n_proc = 10;
c = parcluster('local');
c.NumWorkers = n_proc;
parpool(n_proc);

parfor i=1:length(N)*n_trials

   	fprintf('Test: N = %d, run = %d\n', Nm(i), Run(i));
    
    % generate training data: uniform inputs on [-1,1]
    X = 2*rand(Nm(i), d) - 1;
    Y = func(X);

    % setup structs to save results
    FPCpp  = struct;
    BPDNpp = struct;
    OMPpp  = struct;
    OMPNpp = struct;
    RGApp  = struct;

    % train PC model using LeastSquares 
    FPC  = LeastSquaresPCE(d, order, basis);
    tic; FPC = FPC.fit(X,Y); FPCpp.time = toc;
    [FPCpp.TestE, FPCpp.MeanE, FPCpp.StdE] = FPC.testErr(XTest, YTest, wTest);
    FPCpp.spars = length(FPC.indices);
    
    % train PC model using BPDN 
    BPDN = L1Minimization(d, order, basis);
    tic; BPDN = BPDN.fit(X,Y); BPDNpp.time = toc;
    [BPDNpp.TestE, BPDNpp.MeanE, BPDNpp.StdE] = BPDN.testErr(XTest, YTest, wTest);
    BPDNpp.spars = length(BPDN.indices);

    % train PC model using traditional OMP 
    OMP  = OrthogonalMatchingPursuit(d, order, basis, 'old');
    tic; OMP  = OMP.fit(X,Y); OMPpp.time = toc;
    [OMPpp.TestE, OMPpp.MeanE, OMPpp.StdE] = OMP.testErr(XTest, YTest, wTest);
    OMPpp.spars = length(OMP.indices);

    % train PC model using OMP with optimal basis selection strategy
    OMP  = OrthogonalMatchingPursuit(d, order, basis, 'opt');
    tic; OMP  = OMP.fit(X,Y); OMPNpp.time = toc;
    [OMPNpp.TestE, OMPNpp.MeanE, OMPNpp.StdE] = OMP.testErr(XTest, YTest, wTest);
    OMPNpp.spars = length(OMP.indices);

    % train PC model using Randomized Greedy algorithm
    RGA  = RandomizedGreedy(d, order, basis);
    tic; RGA  = RGA.fit(X,Y); RGApp.time = toc;
    [RGApp.TestE, RGApp.MeanE, RGApp.StdE] = RGA.testErr(XTest, YTest, wTest);
    RGApp.spars = length(RGA.indices);

    % save structs
    FPC_out{i}  = FPCpp;
    BPDN_out{i} = BPDNpp;
    OMP_out{i}  = OMPpp;
    OMPN_out{i} = OMPNpp;
    RGA_out{i}  = RGApp;

end

delete(gcp('nocreate'))

% save matlab file with results
save(['../results/' func_str '_d' num2str(d) '_ord' num2str(order)]);

end
