% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

classdef OrthogonalMatchingPursuit < LeastSquaresPCE
	% OrthogonalMatchingPursuit Class computes coefficients 
	% of PCE expansion given N pairs of samples (X,Y) using
	% the OMP techniques
	%
	% Methods include: fit, train, select_basis, evalCoeffs, 
	% 				   evalResidual, updateQR
	%

	properties

		alg 	 % approach for OMP minimization (opt or old)
		Q 		 % QR decomposition of Psi(:,basis)
		R 		 % QR decomposition of Psi(:,basis)
		PsiNorm  % stored L2 norm of basis functions
		residual % current training residual of approximation

	end

	methods
		function OMP = OrthogonalMatchingPursuit(dim, order, basis, alg, varargin)
			
			% define OMP object
			OMP@LeastSquaresPCE(dim, order, basis, varargin{:});
			OMP.alg = alg;

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function OMP = fit(OMP, X, Y)

			% evaluate basis functions
			Psi = OMP.poly.BasisEval(X);

			% define training function
			solver = @(Xt, Yt, delta) OMP.train(Xt, Yt, delta);

			% run cross-validation to determine \delta parameter
			delta_opt = OMP.cross_validate(solver, Psi, Y);

			% solve for coefficients with optimal delta
			OMP = solver(Psi, Y, delta_opt); 

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function OMP = train(OMP, Psi, Y, delta)

			% save Y, Psi and compute L2 norm of basis evaluations
			OMP.Y        = Y;
			OMP.Psi 	 = Psi;
			OMP.PsiNorm  = sum(Psi.^2,1);
            OMP.residual = Y;
            
            % determine the number of basis elements
			n_basis      = OMP.n_basis();

            % define initial PC dictionary (including constant term)
			OMP.dict 	 = 2:n_basis;
			OMP.indices  = 1;

            % initialize Q, R, and PsiNorm (after updating OMP.indices)
            if strcmp(OMP.alg, 'opt')
                [OMP.Q, OMP.R, OMP.PsiNorm] = OMP.updateQR();
            end

			% initialize L2 norm of residual
            OMP.coeffs   = OMP.evalCoeffs();
            OMP.residual = OMP.evalResidual();
			L2res = norm(OMP.residual, 2);

			% run OMP until dictionary is empty or norm of 
			% training error is less than delta
			while(L2res >= delta && ~isempty(OMP.dict))

				% find next basis_idx and add to index
				basis_idx = OMP.selectBasis();
				OMP.indices = [OMP.indices, basis_idx];
				OMP.dict(OMP.dict == basis_idx) = [];

				% evaluate coefficients
				OMP.coeffs = OMP.evalCoeffs();

				% evaluate training residual
                OMP.residual = OMP.evalResidual();
				L2res = norm(OMP.residual,2);

				% update PsiNorm decomposition
				if strcmp(OMP.alg, 'opt')
					[OMP.Q, OMP.R, OMP.PsiNorm] = OMP.updateQR();
				end

			end

		end	%endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function basis_idx = selectBasis(OMP)

			% extract residual, Psi and PsiNorm
			residual_ = OMP.residual;
			Psi_ 	  = OMP.Psi;
			PsiNorm_  = OMP.PsiNorm;

			% define basis functions to test 
			basis = OMP.dict;

			% extract basis functions and denominators
			Psi_t = Psi_(:,basis);
			PsiNorm_t = PsiNorm_(:,basis);

			% compute correlations of basis vectors with residual
			basis_corr = (residual_'*Psi_t).^2./PsiNorm_t;

			% Identify maximum correlation or assign 1 for the first term
			[~, opt_idx] = max(abs(basis_corr));
			basis_idx = basis(opt_idx);

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function coeffs = evalCoeffs(OMP)
			Psi_ = OMP.Psi(:,OMP.indices);
			coeffs = (Psi_'*Psi_)\(Psi_'*OMP.Y);
		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function residual = evalResidual(OMP)
			residual = OMP.Psi(:,OMP.indices)*OMP.coeffs - OMP.Y;
		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function [Q, R, PsiNorm] = updateQR(OMP)

			% extract Psi
			Psi_ = OMP.Psi(:,OMP.indices);

			% update Q,R
			if isempty(OMP.Q) || isempty(OMP.R)
				[Q,R] = qr(Psi_);
			else
				[Q,R] = qrinsert(OMP.Q, OMP.R, size(Psi_,2), Psi_(:,end));
			end

			% update PsiNorm
            Q_k = Q(:,length(OMP.indices));
			PsiNorm = OMP.PsiNorm - (Q_k'*OMP.Psi).^2;

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
	end %endMethods

end %endClass
