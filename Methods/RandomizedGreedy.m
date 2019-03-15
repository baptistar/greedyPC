% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

classdef RandomizedGreedy < LeastSquaresPCE
	% LeastSquaresPCE Class computes coefficients of PCE
	% expansion given N pairs of samples (X,Y) using
	% approximate least-squares minimization techniques.
	% The implemented methods are: unregularized, OMP, 
	% l1 minimization, and RGA (randomized algorithm).
	%
	% Methods include: fit, post_process, train, select_basis,
	% 				   evalResidual, terminationCriteria, evalCoeffs,
	%				   greedyBacktrack, LeaveOneOut, updateQR
	%

	properties

		Q 		 % QR decomposition of Psi(:,basis)
		R 		 % QR decomposition of Psi(:,basis)
		PsiNorm  % stored L2 norm of basis functions
		residual % current training residual of approximation
		LOO      % track leave one out error for termination criteria

		diffAvg  % parameter for absolute change in mean LOO
		pChange  % parameter for relative change in mean LOO
        RepExp   % parameter for number of repeated RGA experiments
        
	end

	methods
		function RGA = RandomizedGreedy(dim, order, basis, varargin)
			
			% define RGA object
			RGA@LeastSquaresPCE(dim, order, basis, varargin{:});

			% define tolerances for termination criteria
			RGA.diffAvg = 0;
			RGA.pChange = 0.01;
            
            % define parameter for number of RGA experiments
            RGA.RepExp  = 10;

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function RGA = fit(RGA, X, Y)

			% evaluate basis functions
			Psi = RGA.poly.BasisEval(X);

			% repeat RGA training for RepExp times
            RGA_PCE = cell(RGA.RepExp,1);
            for i=1:RGA.RepExp
                RGA_PCE{i} = RGA.train(Psi, Y);
            end
            
            % extract LOO errors
            RGA_LOO = zeros(RGA.RepExp,1);
            for i=1:RGA.RepExp
                RGA_LOO(i) = min(RGA_PCE{i}.LOO);
            end
            
            % find expansion with minimum LOO error
            [~,opt_idx] = min(RGA_LOO);
            RGA = RGA_PCE{opt_idx};

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function RGA = train(RGA, Psi, Y)

			% save Y, Psi and compute L2 norm of basis evaluations
			RGA.Y        = Y;
			RGA.Psi 	 = Psi;
			RGA.PsiNorm  = sum(Psi.^2,1);
            RGA.residual = Y;

            % determine the number of basis elements
			n_basis = RGA.n_basis();

            % define PC dictionary (including constant term)
			RGA.dict = 2:n_basis;
			RGA.indices = 1;
				
            % initialize Q, R, and PsiNorm (after updating RGA.indices)
            [RGA.Q, RGA.R, RGA.PsiNorm] = RGA.updateQR();

            % update residual and LOO error (after updating RGA.indices)
            RGA.residual = RGA.evalResidual();
            RGA.LOO      = RGA.LeaveOneOut();

			% initialize termination criteria
			term_crit = 0;

			while(term_crit == 0 && ~isempty(RGA.dict))

				% select next basis function
				basis_idx = RGA.select_basis();
				RGA.indices = [RGA.indices, basis_idx];
				RGA.dict(RGA.dict == basis_idx) = [];

				% update Q, R, and PsiNorm
				[RGA.Q, RGA.R, RGA.PsiNorm] = RGA.updateQR();

				% evaluate residual and Leave-One-Out
				RGA.residual = RGA.evalResidual();
				RGA.LOO = [RGA.LOO, RGA.LeaveOneOut()];

				% evaluate termination criteria
				term_crit = RGA.terminationCriteria();
			end

			% backtrack solution to minimum value of stopping criteria
			RGA = RGA.greedyBacktrack();
			
			% evaluate coefficients and post_process solution
			RGA.coeffs = RGA.evalCoeffs();

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function basis_idx = select_basis(RGA)

			% extract residual, Psi and PsiNorm
			residual_ = RGA.residual;
			Psi_ 	  = RGA.Psi;
			PsiNorm_  = RGA.PsiNorm;

			% define basis functions to test 
			n_basis = min(60,length(RGA.dict));
			[basis, ~] = datasample(RGA.dict, n_basis, 'Replace', false);

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
		function [residual] = evalResidual(RGA)
			% update residual using QR decomposition
            Q_k_ = RGA.Q(:,length(RGA.indices));
			residual = RGA.residual - Q_k_'*RGA.Y*Q_k_;
		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function term_crit = terminationCriteria(RGA)

			% compute moving averages of LOO starting at iteration 10
			if length(RGA.indices) >= 10

				% compute mean over two time frames
				mu_1 = mean(RGA.LOO(end-4:end));
				mu_2 = mean(RGA.LOO(end-9:end-5));
			
				% compute decrease and percent change in LOO error
				diff_avg = mu_1 - mu_2;
				comp_avg = abs((mu_1 - mu_2)/mu_2);

			else
				diff_avg = nan;
				comp_avg = nan;

			end

			% evaluate termination crteria
			term_crit = (diff_avg > RGA.diffAvg) || ...
						(comp_avg < RGA.pChange) || ...
						(length(RGA.indices) == size(RGA.Psi,1));

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function coeffs = evalCoeffs(RGA)
            Q_ = RGA.Q(:,1:length(RGA.indices));
            R_ = RGA.R(1:length(RGA.indices),:);
			coeffs = R_\(Q_'*RGA.Y);
		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function RGA = greedyBacktrack(RGA)
        % update properties of PC approximation (Q, R, PsiNorm, residual)
        % to setting at minimum of LOO error 
            
			% find minimum of stopping criteria (LOO)
			[~, min_idx] = min(RGA.LOO);

			% backtrack removing basis functions
			for i=length(RGA.indices):-1:(min_idx+1)

                % update PsiNorm and residual
                Q_k_ = RGA.Q(:,length(RGA.indices));
				RGA.PsiNorm  = RGA.PsiNorm + (Q_k_'*RGA.Psi).^2;
				RGA.residual = RGA.residual + Q_k_'*RGA.Y*Q_k_;

				% udate QR decomposition (by removing column)
				[RGA.Q, RGA.R] = qrdelete(RGA.Q, RGA.R, i);
                
				% update dictionary and indices
				RGA.dict = sort([RGA.dict, RGA.indices(end)]);
				RGA.indices(end) = [];

			end

		end %endFunction
        %------------------------------------------------------------------
		%------------------------------------------------------------------
		function LOO = LeaveOneOut(RGA)

            % extract residual and Q_k
            residual_ = RGA.residual;
            Q_ = RGA.Q(:,1:length(RGA.indices));
            
			% compute gradient with entrywise min for stability
			res_gradient = 1 - sum(Q_.^2,2);
			res_gradient = max(abs(res_gradient), 1e-12);

			% compute LOO error
			N = size(RGA.Psi,1);
			LOO = 1/N*sum((residual_./res_gradient).^2);

		end %endFunction
   		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function [Q, R, PsiNorm] = updateQR(RGA)
        % update Q, R and PsiNorm after adding a new index to Psi 

			% extract Psi
			Psi_ = RGA.Psi(:,RGA.indices);

			% update Q,R
			if isempty(RGA.Q) || isempty(RGA.R)
				[Q,R] = qr(Psi_);
			else
				[Q,R] = qrinsert(RGA.Q, RGA.R, size(Psi_,2), Psi_(:,end));
			end

			% update PsiNorm
            Q_k_ = Q(:,length(RGA.indices));
			PsiNorm = RGA.PsiNorm - (Q_k_'*RGA.Psi).^2;

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
	end %endMethods

end %endClass
