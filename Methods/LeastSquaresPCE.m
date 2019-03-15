% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

classdef LeastSquaresPCE
	% LeastSquaresPCE Class computes coefficients of PCE
	% expansion given N pairs of samples (X,Y) using
	% approximate least-squares minimization techniques.
	% 
	% The implemented methods defined as subclasses are:
	%  - unregularized L2 minimization (parent)
	%  - OMP & OMP with optimal basis selection
	%  - L1 minimization (BPDN)
	%  - Randomized greedy (RGA)
	%
	% Methods include: cross_validate, validation_error, fit, testErr, n_basis
	%

	properties

		dim    		% input dimension
		order 		% order of expansion
		poly        % multivariate polynomial object

		dict        % list of available basis functions
		indices     % list of selected basis functions
		coeffs      % basis coefficients

		Psi         % evaluations of basis functions
		Y 			% evaluations of model output

    end

	methods
		function PC = LeastSquaresPCE(dim, order, basis, varargin)
			
			% Define PC object
			p = ImprovedInputParser;
			addRequired(p,'dim');
			addRequired(p,'order');
			parse(p,dim,order,varargin{:});
			PC = passMatchedArgsToProperties(p, PC);

			% save multivariate PC object
			PC.poly = MultivariatePC(dim, order, basis);

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function nPC = n_basis(PC)
			nPC = PC.poly.ncoeff();
		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function PC = fit(PC,X,Y)
			
			% evaluate basis functions
			Psi_ = PC.poly.BasisEval(X);

			% save Y, Psi
			PC.Y    = Y;
			PC.Psi 	= Psi_;

			% evaluate coefficients using normal equations
			PC.coeffs = (Psi_'*Psi_)\(Psi_'*Y);

			% define indices and empty dictionary
			PC.indices = 1:PC.n_basis();
			PC.dict    = [];

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function delta = cross_validate(PC, PCEsolver, X, Y)
		
			% define number of Cross-Validation folds
			n_folds = 4;

			% residual error evaluations
			delta_cv = logspace(-5,1,100);

			% determine total number of samples
			N = size(X,1);

            % define matrix to store test error
            test_error = nan(n_folds, length(delta_cv));
            
			for i=1:n_folds

				% Determine the number of terms in each set
				pts_pfold  = floor(N/n_folds);
				test_pts   = (i-1)*pts_pfold+1:i*pts_pfold;
				train_pts  = setdiff(1:N, test_pts);

				% separate dataset into training and test sets
				xTrain = X(train_pts,:);
				yTrain = Y(train_pts,:);

				xTest  = X(test_pts,:);
				yTest  = Y(test_pts,:);

				% Train model and determine the training-set error
				for j=1:length(delta_cv)
					PCE = PCEsolver(xTrain, yTrain, delta_cv(j));
					test_error(i,j) = PCE.validation_error(xTest, yTest);
				end

			end

			% compute average test error
			mean_error = nanmean(test_error,1);

			% Find minimum validation error
			[~, idx_min] = min(mean_error);
			delta = delta_cv(idx_min(1));

			% Rescale optimal delta based on number of samples
			% consider minimum of delta in case of multiple values
			delta = sqrt(N/length(train_pts))*min(delta);

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function err = validation_error(PC, Psi, Y)
			% evaluate L2 norm of prediction error at (X,Y)
			err = norm(Psi(:,PC.indices)*PC.coeffs - Y, 2);
		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function [TestE, MeanE, StdE] = testErr(PC, XTest, YTest, wTest)

			% evaluate basis functions at XTest
			Psi_ = PC.poly.BasisEval(XTest);
			Yhat = Psi_(:,PC.indices)*PC.coeffs;

			% compute relative weighted test error
			L2Err = (YTest - Yhat).^2;
			TestE = (wTest*L2Err)/(wTest*(YTest.^2));

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
			MeanE = norm(meanPC - meanTest)/norm(meanTest);
			StdE  = norm(stdPC - stdTest)/norm(stdTest);

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
	end %endMethods

end %endClass
