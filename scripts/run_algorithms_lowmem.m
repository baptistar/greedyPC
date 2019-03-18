% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

function run_algorithms_lowmem(test)
% RUN_ALGORITHMS_LOWMEM: Runs the linear regression algorithms
% on the test case specified in the test struct without any
% sparse grid post_processing. The results are saved in the 
% results/test_name file.

% extract parameters
d          = test.d;
order      = test.order;
basis      = test.basis;
func       = test.func;
func_str   = test.func_str;
N          = test.N;
n_trials   = test.n_trials;

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
    FPCpp.spars = length(FPC.indices);
    FPCpp.PCE = FPC;    

    % train PC model using BPDN 
    BPDN = L1Minimization(d, order, basis);
    tic; BPDN = BPDN.fit(X,Y); BPDNpp.time = toc;
    BPDNpp.spars = length(BPDN.indices);
    BPDNpp.PCE = BPDN;    

    % train PC model using traditional OMP 
    OMP  = OrthogonalMatchingPursuit(d, order, basis, 'old');
    tic; OMP  = OMP.fit(X,Y); OMPpp.time = toc;
    OMPpp.spars = length(OMP.indices);
    OMPpp.PCE = OMP;    

    % train PC model using OMP with optimal basis selection strategy
    OMP  = OrthogonalMatchingPursuit(d, order, basis, 'opt');
    tic; OMP  = OMP.fit(X,Y); OMPNpp.time = toc;
    OMPNpp.spars = length(OMP.indices);
    OMPNpp.PCE = OMP;

    % train PC model using Randomized Greedy algorithm
    RGA  = RandomizedGreedy(d, order, basis);
    tic; RGA  = RGA.fit(X,Y); RGApp.time = toc;
    RGApp.spars = length(RGA.indices);
    RGApp.PCE = RGA;

    % save structs
    FPC_out{i}  = FPCpp;
    BPDN_out{i} = BPDNpp;
    OMP_out{i}  = OMPpp;
    OMPN_out{i} = OMPNpp;
    RGA_out{i}  = RGApp;

end

delete(gcp('nocreate'))

% save matlab file with results
save(['../results/' func_str '_d' num2str(d) '_ord' num2str(order) '_lowmemory']);

% run post-processing
post_process_lowmem(test)

end
