% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

classdef L1Minimization < LeastSquaresPCE
	% L1Minimization Class computes coefficients 
	% of PCE expansion given N pairs of samples (X,Y) using
	% the BPDN L1 minimization algorithm
	%
	% Methods include: fit, train, eval_residual
	%
	
	properties

	end

	methods
		function L1 = L1Minimization(dim, order, basis, varargin)
			
			% define LS PC object
			L1@LeastSquaresPCE(dim, order, basis, varargin{:});

			% check if spg_bpdn is in the path
			if exist('spg_bpdn','file') ~= 2
				error('Please add SPGL1 solver to path - see README');
			end

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function L1 = fit(L1, X, Y)

			% evaluate basis functions
			Psi = L1.poly.BasisEval(X);

			% save Y, and Psi
			L1.Y    = Y;
			L1.Psi 	= Psi;

			% define training function
			solver = @(Xt, Yt, delta) L1.train(Xt, Yt, delta);

			% run cross-validation to determine \delta parameter
			delta_opt = L1.cross_validate(solver, Psi, Y);

			% solve for coefficients with optimal delta
			L1 = solver(Psi, Y, delta_opt); 

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function L1 = train(L1, Psi, Y, delta)

			% extract normalization of basis to define preconditoner
			precond = sqrt(sum(Psi.^2,1));
			W = diag(1./precond);

			% define A, b matrices for BPDN solver
			A = Psi*W;
			b = Y;

			% Solve for the solution
			options = spgSetParms('verbosity',0,'weights',diag(W));
			[soln,~,~,info] = spg_bpdn(A, b, delta, options);

			% If not converged, print error message
			if (info.stat > 5)
				warning('BPDN solver did not converge.')
			end

			% save solution in L1 structure
			soln = W*soln;
			L1.coeffs  = nonzeros(soln);
			L1.indices = find(soln);
			L1.dict    = setdiff(1:L1.n_basis(), L1.indices);

		end	
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function residual = eval_residual(L1)
			residual = L1.Psi(:,L1.indices)*L1.coeffs - L1.Y;
		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
	end %endMethods

end %endClass
