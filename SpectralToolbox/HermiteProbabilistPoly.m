% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

classdef HermiteProbabilistPoly
	% HermiteProbabilistPoly defines Probabilist Hermite polynomials
	% using the rescaling of the physicist Hermite polynomials
	%
	% Un-normalized:
	%  He_0(x) = 1
	%  He_1(x) = x
	%  He_2(x) = x^2 - 1
	%
	% Normalized:
	%  \tilde{He}_n(x) = He_n(x)/\sqrt(\sqrt(2*pi)*n!)
	%
	% Methods include: Evaluate, GradEvaluate

	properties

		HPhy  % Physicist Hermite polynomial object

	end

	methods
		function HP = HermiteProbabilistPoly(varargin)
			
			% Define HP object
			p = ImprovedInputParser;
			parse(p, varargin{:});
			HP = passMatchedArgsToProperties(p, HP);

			% Set Physicist polynomial property
			HP.HPhy = HermitePhysicistPoly();
		
		end %endFunction
		%---------------------------------------------------------
		%--------------------------------------------------------- 
		function p = Evaluate(HP, x, N, norm)
		% Evaluate the N-th order Hermite probabilists' polynomial
		% using the relation He_n(x) = 2^(-n/2)*H_n(x/sqrt(2))
		% accounting for the correct scaling factors
			
			p = 1/sqrt(2^N)*HP.HPhy.Evaluate(x/sqrt(2),N,false);
			if norm == true
				p = p/sqrt(sqrt(2*pi) * gamma(N+1));
			end

		end %endFunction
		%---------------------------------------------------------
		%--------------------------------------------------------- 
		function dp = GradEvaluate(HP, x, N, k, norm)
		% Evaluate the k-th derivative of the N-th order Hermite 
		% probabilists' polynomial at m=size(x,1) points using the
		% relation He_n(x) = 2^(-n/2)*H_n(x/sqrt(2)) and accounting
		% for the correct scaling factors

			dp = 1/sqrt(2^N)*HP.HPhy.GradEvaluate(x/sqrt(2),N,k,false);
			dp = dp*(1/sqrt(2))^k;
			if norm == true
				dp = dp/sqrt(sqrt(2*pi) * gamma(N+1));
			end

		end %endFunction
		%---------------------------------------------------------
		%--------------------------------------------------------- 
	end %endMethods

end %endClass
