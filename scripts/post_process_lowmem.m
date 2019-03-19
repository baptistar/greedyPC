% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

function post_process_lowmem(test)
% POST_PROCESS_LOWMEM: Runs the post processing using the 
% the sparse grid model evaluations for the FPC, BPDN, OMP
% OMPN and RGA algorithms after running the script
% run_algorithms_lowmem.m

% extract parameters
d          = test.d;
order      = test.order;
grid_level = test.grid_level;
func       = test.func;
func_str   = test.func_str;
N          = test.N;
n_trials   = test.n_trials;

% load data from file
file_name = ['../results/' func_str '_d' num2str(d) '_ord' num2str(order) '_lowmemory'];
load(file_name,'FPC_out','BPDN_out','OMP_out','OMPN_out','RGA_out');
models = {FPC_out, BPDN_out, OMP_out, OMPN_out, RGA_out};

% call sparse grid to determine the test samples
[XTest, wTest] = sparse_grid(d, grid_level, 'unif');
YTest = func(XTest);

% evaluate all basis functions at XTest
Psi_ = FPC_out{1,1}.PCE.poly.BasisEval(XTest);

for k=1:length(models)

	% extract model
	model_out = models{k};

	for i=1:length(N)
		for j=1:n_trials

			fprintf('Test: N = %d, run = %d\n', N(i), j);

			% extract PC expansion
			PC = model_out{i,j}.PCE;
		
			% evaluate expansion at XTest
			coeff = zeros(size(Psi_,2),1);
			coeff(PC.indices) = PC.coeffs;
			Yhat = Psi_*coeff;

			% compute relative weighted test error
			L2Err = (YTest - Yhat).^2;
			model_out{i,j}.TestE = (wTest*L2Err)/(wTest*(YTest.^2));

			% evaluate mean and standard deviation of yTest
			meanTest = wTest*YTest;
			stdTest  = sqrt(wTest*YTest.^2 - meanTest^2);

			% evaluate mean and standard deviation of PC expansion
			if PC.indices(1) == 1
				meanPC = PC.coeffs(1);
			else
				meanPC = mean(Yhat);
			end
			stdPC = sqrt(sum((PC.coeffs).^2) - meanPC^2);

			% compute relative mean and standard deviation error
			model_out{i,j}.MeanE = norm(meanPC - meanTest)/norm(meanTest);
			model_out{i,j}.StdE  = norm(stdPC - stdTest)/norm(stdTest);

			model_out{i,j}.PCE = [];

		end
	end
	
	% save model_out
	models{k} = model_out;

end

% extract models
FPC_out  = models{1};
BPDN_out = models{2};
OMP_out  = models{3};
OMPN_out = models{4};
RGA_out  = models{5};

% save matlab file with results
save(['../results/' func_str '_d' num2str(d) '_ord' num2str(order)]);

end
