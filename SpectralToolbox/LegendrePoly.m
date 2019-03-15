% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

classdef LegendrePoly
	% LegendrePoly defines Legendre polynomials on [-1,1]
	%
	% Un-normalized:
	%  H_0(x) = 1
	%  H_1(x) = x
	%  H_2(x) = 0.5*(3x^2 - 1)
	%
	% Normalized:
	%  \tilde{P}_n(x) = H_n(x)/\sqrt(2 / (2*n+1))
	%
	% Methods include: Evaluate, MonomialCoeffs, GradEvaluate

	properties

	end %endProperties

	methods
		function LP = LegendrePoly(varargin)
			
			% Define LP object
			p = ImprovedInputParser;
			parse(p, varargin{:});
			LP = passMatchedArgsToProperties(p, LP);

		end %endFunction
		%---------------------------------------------------------
		%--------------------------------------------------------- 
		function c = MonomialCoeffs(LP, N, norm)
		% Return the monomial polynomial coefficients.
		% Normalized coefficients divide by \sqrt(1/(2*N+1))
			c = LP.LegCoeffs(N);
			if norm == true
				c = c/sqrt(1/(2*N+1));
			end
		end %endFunction
		%---------------------------------------------------------
		%--------------------------------------------------------- 
		function c = LegCoeffs(LP, N)
		% Compute the coefficients of the Legendre polynomials
		% using a recursive algorithm. Input: N = polynomial order
		%
		% Example: [1], P(0,x) = 1
		%		   [0, 1], P(1,x) = x
		% 		   [-0.5, 0, 1.5], P(2,x) = 1.5x^2 - 0.5; 
		%
		% Reference: adapted from John Burkardt (2004)

			% N<0, N=0, and N=1 cases
			if N < 0
				c = [];
				return
			elseif N == 0
				c = 1;
				return
			elseif N == 1
				c = [0 1];
				return
			end

			% compute N>=2 by recursion
			c = [1, 0; 0, 1];
			for i=2:N
				c(i+1,1:i-1) =                (   - i + 1 ) * c(i-1,1:i-1) / (i);
				c(i+1,2:i+1) = c(i+1,2:i+1) + ( i + i - 1 ) * c(i  ,1:i  ) / (i);
			end
			c = c(end,:);
		
		end %endFunction
		%---------------------------------------------------------
		%--------------------------------------------------------- 
		function p = Evaluate(LP, x, N, norm)
		% Evaluate the N-th order Legendre polynomial at m=size(x,1) 
		% points using the recursion relation P_0(x) = 1 and
		% nP_{n}(x) = (2n-1)*x*P_{n-1}(x) - (n-1)P_{n-2}(x)

			% set initial condition P_0(x) for recursion
			m = size(x,1);
			p = ones(m,1);

			% run recursion by tracking P_{n-1}, P_{n}, P_{n+1}
			if N > 0

				% initialize p_jm1, p_jm2
				p_jm1 = p;
				p_jm2 = zeros(m,1);

				% evaluate H_{j} by recursion and
				% saving H_{j-1}(x) and H_{j}(x)
				for j=1:N
					p = (2*j-1)/j*x.*p_jm1 - (j-1)/j*p_jm2;
					p_jm2 = p_jm1;
					p_jm1 = p;
				end
			end

			% normalize polynomials
			if norm == true
				p = p/sqrt(1/(2*N+1));
			end

		end %endFunction
		%---------------------------------------------------------
		%--------------------------------------------------------- 
		function dp = GradEvaluate(LP, x, N, k, norm)
		% Evaluate the k-th derivative of the N-th order Legendre 
		% polynomials at m=size(x,1) points using the
		% relation (x^2-1)/N P_{n}'(x) = xP_{n}(x) - P_{n-1}(x)

			if k == 0
				dp = LP.Evaluate(x, N, norm);
			elseif k == 1
				dp = x.*LP.Evaluate(x, N, false) - ...
					LP.Evaluate(x, N-1, false);
				dp = N./(x.^2 - 1).*dp;
				if norm == true
					dp = dp/sqrt(1/(2*N+1));
				end
			else
				error('Derivative is not implemented for k > 1.')
			end

		end %endFunction
		%---------------------------------------------------------
		%--------------------------------------------------------- 
	end %endMethods

end %endClass
