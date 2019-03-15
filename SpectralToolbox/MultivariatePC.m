% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

classdef MultivariatePC
	% MultivariatePC defines a class of multivariate
	% Hermite or Legendre polynomial basis functions: 
	% 	H(x) = \sum_{\alpha} \Psi_{\alpha}*c_{\alpha}
	%
	% Methods: BasisEval, dxBasisEval, Evaluate, Grad_a, Grad_x
	%
	% Author: Ricardo Baptista
	% Date:   November 2018

	properties

		poly            % object for univariate base polynomial
		dim    			% input dimension
		order  			% polynomial order
		norm            % normalization property of monomials
		coeff  		    % total order expansion coefficients

		basis_precomp   % precomputed basis evaluations
		dxbasis_precomp % precomputed derivatives of basis evaluations

	end

	methods
		function mPC = MultivariatePC(dim,order,basis,varargin)

			% Define TM object
			p = ImprovedInputParser;
			addRequired(p,'dim');
			addRequired(p,'order');
			parse(p,dim,order,varargin{:});
			mPC = passMatchedArgsToProperties(p, mPC);

			% Define hermite base polynomial
			if strcmp(basis,'Legendre')
				mPC.poly = LegendrePoly();
			elseif strcmp(basis, 'Hermite')
				mPC.poly = HermiteProbabilistPoly();
			else
				error('Basis is not implemented')
			end
			
			% set normalization property
			mPC.norm = true;

			% Initialize precomputed
			mPC.basis_precomp   = [];
			mPC.dxbasis_precomp = [];

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function d = ncoeff(mPC)
			d = nchoosek(mPC.dim + mPC.order, mPC.order);
		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function mIndices = totalDeg(mPC, d, order)
		% totalDeg: Function computes indices of total degree PC 
		% expansion. Output is a cell with indices for each order
		%
		% Reference: MATLAB polynomial chaos toolbox

			% declare cell to store indices for each order
			mIndices    = cell(order+1,1); % multi-index
			mIndices{1} = zeros(1,d);      % multi-index for length 0

			if d == 1
				for q=1:order
					mIndices{q+1} = q;
				end
			else
				for q = 1:order
					s  = nchoosek(1:d+q-1,d-1);
					s1 = zeros(size(s,1),1);
					s2 = (d+q)+s1;
					mIndices{q+1} = flipud(diff([s1 s s2],1,2))-1; % -1 due to MATLAB indexing
					if sum(mIndices{q+1},2) ~= q*ones(nchoosek(d+q-1,d-1),1)
						error('The sum of each row has to be equal to q-th order');
					end
				end
			end

			% convert orders to matrix
			mIndices = cell2mat(mIndices);

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function Psi = BasisEval(mPC, X)
			
			% check if data is empty and return constant in this case
			m = size(X,2);
            if m == 0
				Psi = mPC.poly.Evaluate(ones(m,1), 0, mPC.norm);
            else
                Psi = mPC.PolyVandermonde(X, 0, []);
            end

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function dxPsi = dxBasisEval(mPC, X)
			dxPsi = cell(mPC.d,1);
			for i=1:d
				dxPsi{i} = mPC.PolyVandermonde(X, 1, i);
			end
			dxPsi = cell2mat(dxPsi);
		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function Psi = PolyVandermonde(mPC, X, k, grad_dim)

			% determine the number of samples and parameters
			[N, d] = size(X);

			% find multi-indices for total degree PC expansion
			mIndices = mPC.totalDeg(d, mPC.order);
			n_basis  = size(mIndices,1);

			% declare matrix to store polynomial evaluations
			Psi = zeros(N, n_basis);

			% evaluate all basis functions
			for i=1:n_basis

				% extract row from basis_prod
				poly_ind = mIndices(i,:);

				% extract the value of the basis functions for all data values
				poly_all = zeros(N, d);
				for j=1:d
					if (k ~= 0) && (j==grad_dim)
						poly_all(:,j) = mPC.poly.GradEvaluate(X(:,j), poly_ind(j), k, mPC.norm);
					else
						poly_all(:,j) = mPC.poly.Evaluate(X(:,j), poly_ind(j), mPC.norm);
					end
				end

				% compute product of 1D basis functions and assign to PC_basis_value
				Psi(:,i) = prod(poly_all,2);

			end

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function S = Evaluate(mPC, X)

			% check dimension of input samples
			if size(X,2) ~= mPC.dim
				error('Inputs and poly dimension mismatch')
            end
            
			% evaluate S(x) = Psi(x)*coeff
			if isempty(mPC.basis_precomp)
				mPC.basis_precomp = mPC.BasisEval(X);
            end
            S = mPC.basis_precomp*mPC.coeff;

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function daS = Grad_a(mPC, X)

			% check dimension of input samples
			if size(X,2) ~= mPC.dim
				error('Inputs and poly dimension mismatch')
			end

			% evaluate \nabla_coeff S(x) = Psi(x)
			if isempty(mPC.basis_precomp)
				mPC.basis_precomp = mPC.BasisEval(X);
			end
			daS = mPC.basis_precomp;

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function dxS = Grad_x(mPC, X)

			% check dimension of input samples
			if size(X,2) ~= mPC.dim
				error('Inputs and poly dimension mismatch')
			end

			% evaluate (\nabla_x Psi(x))mPC.coeff,1,1,mPC.ncoeff
            if isempty(mPC.dxbasis_precomp)
				mPC.dxbasis_precomp = mPC.dxBasisEval(X);
            end

			% evaluate \nabla_x S(x) = (\nabla_x Psi(x))*coeff with reshape
            coeff_rep = repmat(reshape(mPC.coeff,1,1,mPC.ncoeff),size(X,1),size(X,2),1);
			dxS = sum(mPC.dxbasis_precomp.*coeff_rep,3);

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function d2aS = Hess_a(mPC, X)

			% find number of input samples
			N = size(X,1);

			% evaluate \nabla^2_coeff S(x) = 0
			d2aS = zeros(N, mPC.ncoeff, mPC.ncoeff);

		end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
	end %endMethods

end %endClass
