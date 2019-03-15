% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

classdef HermitePhysicistPoly
	% HermitePhysicistPoly defines Physicist Hermite polynomials 
	%
	% Un-normalized:
	%  H_0(x) = 1
	%  H_1(x) = 2*x
	%  H_2(x) = 4x^2 - 2
	%
	% Normalized:
	%  \tilde{H}_n(x) = H_n(x)/\sqrt(n!*2^n*\sqrt(pi))
	%
	% Methods include: Evaluate, GradEvaluate

	properties

	end %endProperties

	methods
		function HP = HermitePhysicistPoly(varargin)
			
			% Define HP object
			p = ImprovedInputParser;
			parse(p, varargin{:});
			HP = passMatchedArgsToProperties(p, HP);

		end %endFunction
		%---------------------------------------------------------
		%--------------------------------------------------------- 
		function c = MonomialCoeffs(HP, N, norm)
		% Return the monomial polynomial coefficients
		% normalized coefficients divide by \sqrt(\sqrt(pi) * 2^N * N!)
			c = HP.HermiteCoeffs(N);
			if norm == true
				c = c/sqrt(sqrt(pi) * 2^N * factorial(N));
			end
		end %endFunction
		%---------------------------------------------------------
		%--------------------------------------------------------- 
		function cf = HermiteCoeffs(HP, N)
		% Compute the coefficients of the Hermite polynomials
		% using a recursive algorithm. Input: N = polynomial order
		%
		% Example: [1], H(0,x) = 1
		%		   [0, 2], P(1,x) = 2*x
		% 		   [-2, 0, 4], P(2,x) = 4x^2 - 2; 
		%
		% Reference: adapted from John Burkardt (2010)

			% N<0, N=0, and N=1 cases
			if N < 0
				c = [];
				return
			elseif N == 0
				c = 1;
				return
			elseif if N == 1
				c = [0 2];
				return
			end

			% compute N>=2 by recursion
			c = [1, 0; 0, 2];
			for i=2:N
				c(i+1,1)     =  -2.0 * ( i - 1 ) * c(i-1,1);
				c(i+1,2:i-1) =   2.0             * c(i  ,1:i-2)...
								-2.0 * ( i - 1 ) * c(i-1,2:i-1);
				c(i+1,  i  ) =   2.0             * c(i  ,  i-1);
				c(i+1,  i+1) =   2.0             * c(i  ,  i  );
			end
			c = c(end,:);
		
		end % endFunction
		%---------------------------------------------------------
		%--------------------------------------------------------- 
		function p = Evaluate(HP, x, N, norm)
		% Evaluate the N-th order Hermite physicist polynomial
		% at m=size(x,1) points using the recursion relation 
		% H_0(x) = 1 and H_{n}(x) = 2xH_{n-1}(x) - H_{n-1}'(x)
		% where H_{n-1}'(x) = 2(n-1)H_{n-2}(x)
			
			% set initial condition H_0(x) for recursion
			m = size(x,1);
			p = ones(m,1);

			% run recursion by tracking H_{n-1}, H_{n}, H_{n+1}
			if N > 0

				% initialize p_jm1, p_jm2
				p_jm1 = p;
				p_jm2 = zeros(m,1);

				% evaluate H_{j} by recursion and
				% saving H_{j-1}(x) and H_{j}(x)
				for j=1:N
					p = 2*x.*p_jm1 - 2*(j-1)*p_jm2;
					p_jm2 = p_jm1;
					p_jm1 = p;
				end
			end

			% normalize polynomials
			if norm == true
				p = p/sqrt(gamma(N+1) * 2.^N * sqrt(pi));
			end

		end %endFunction
		%---------------------------------------------------------
		%--------------------------------------------------------- 
		function dp = GradEvaluate(HP, x, N, k, norm)
		% Evaluate the k-th derivative of the N-th order Hermite 
		% physicist polynomial at m=size(x,1) points using the
		% relation H_{n}^(k)(x) = 2^{k} n!/(n-k)! H_{n-k}(x)
			
			if N >= k
				fact = 2^k*exp(gammaln(N+1) - gammaln(N-k+1));
				dp = fact*HP.Evaluate(x, N-k, norm);
				if norm == true
					dp = dp/sqrt(fact);
				end
			else
				dp = zeros(size(x));
			end
		
		end %endFunction
		%---------------------------------------------------------
		%--------------------------------------------------------- 
	end %endMethods

end %endClass
