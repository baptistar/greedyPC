% Author: Ricardo Baptista and Prasanth Nair
% Date:   March 2019
%
% See LICENSE.md for copyright information
%

function [global_M, global_K, global_F, KL, p, e, t, int_p] = ...
            diff_2D_stochastic_forcing(n_par, sigma, lc, dim_lim)
% DIFF_2D_STOCHASTIC_FORCING: Function assembles the 2D discretized
% weak form of the parametrized diffusion equation based on a 
% Karhunen-Loeve expansion with a Gaussian correlation function 
% for the parameters
%
% Inputs: n_par   - number of parameters in the problem
%         mean_d  - mean of parameters
%         sigma_d - variance of parameters
%         lc      - correlation length of parameters
%         dim_lim - dimensions of spatial coordinates

%% Geometry and Boundary Conditions

% Define rectangular geometry
g = decsg([3 4 0 1 1 0 0 0 1 1]');

% Define the mesh for the geometry
[p,e,t] = initmesh(g);
[p,e,t] = refinemesh(g,p,e,t);
fprintf('Created a mesh with %i nodes.\n', length(p));

% Get the indexes of interior/boundary FE nodes
[ind_interior, ~] = getindexesMesh(p, dim_lim);

% Set values of interior nodes
int_p = p(:,ind_interior);

%% Initial Solve of Deterministic PDE

% Coefficients of elliptic PDE problem -div(c*grad(u))+a*u=f
c = 1;
a = 0;

% Define the source term as a spatial function
sf = '-1';

% Extract the Mass and Stiffness Matrices
[K, ~, F] = assema(p,t,c,a,sf); % M_other is NULL
[~, M, ~] = assema(p,t,c,1,sf); % M is a SPD matrix

%% KL Expansion

% Declare matrix for correlation evaluations
szmesh=size(p,2);
C = zeros(szmesh, szmesh);

% Evaluate the correlation function
for i=1:szmesh
  for j=1:szmesh
     C(i,j) = exp(-(p(1,i)-p(1,j))^2/(2*lc(1)^2)-(p(2,i)-p(2,j))^2/(2*lc(2)^2));
  end
end
C = sigma^2*C;

% Find the eigenvalues of the KL expansion
W = M*C*M;
W = 0.5*(W+W'); % Symmetrization (otherwise W is not perfectly symmetric) 
[VV,DD] = eig(W,full(M)); 

% Check that the eigenvectors are orthonormalized wrt the mass matrix:
% VV'*M*VV

% Extract the eigenvalues and sort by increasing order
diagDD = diag(DD);
eigenvalKL = sort(diagDD);

% Extract the spatial modes of the KL expansion
ksort = diag(eigenvalKL);
spModesKL = VV(:,szmesh-n_par+1:szmesh)*...
            sqrt(ksort(szmesh-n_par+1:szmesh,szmesh-n_par+1:szmesh));
     
% Extract the KL expansion
KL = spModesKL(ind_interior,:);
        
% Extract mesh parameters
nbTriangles = size(t,2);
A = ones(3,3);
RFcentroids = zeros(nbTriangles,n_par);

% Interpolate the Basis Functions on the Triangular Mesh
for m=1:n_par
   RFm=spModesKL(:,m);
   for k=1:nbTriangles
      ii=t(1:3,k);
      alpha=RFm(ii);
      coordCentroid(1)=mean(p(1,ii)); 
      coordCentroid(2)=mean(p(2,ii));
      A(:,2)=p(1,ii)'; 
      A(:,3)=p(2,ii)';
      beta=A\alpha;
      RFcentroids(k,m)=beta(1)+beta(2)*coordCentroid(1)+beta(3)*coordCentroid(2);
   end
end

%% Final Solution of PDE

% Extract stiffness matrices for KL expansion
KKL = cell(n_par,1);
for m=1:n_par
    cm = RFcentroids(:,n_par-m+1)';
    [KKL{m}, ~, ~] = assema(p, t, cm, 0, zeros(nbTriangles,1)');
end

% Extract right hand side vectors for KL expansion
FKL = cell(n_par,1);
for m=1:n_par
    cm = RFcentroids(:,n_par-m+1)';
    [~, ~, FKL{m}] = assema(p, t, 0, 0, cm);
end

% Extract mass/stiffness matrix elements for interior nodes
M_int = M(ind_interior,ind_interior); % Mass Matrix
K_int = K(ind_interior,ind_interior); % (Main) Stiffness Matrix
F_int = F(ind_interior);

% Extract KL stiffness matrix elements for interior nodes
KKL_int = cell(n_par, 1);
FKL_int = cell(n_par, 1);
for m = 1:n_par
    KKL_int{m} = KKL{m}(ind_interior,ind_interior);
    FKL_int{m} = FKL{m}(ind_interior);
end

% Assemble terms for function call
global_M = M_int;
global_K = [{K_int}; KKL_int];
global_F = [{F_int}; FKL_int];

end

function [ind_interior,ind_boundary] = getindexesMesh(p, dim_lim)

  xmin=dim_lim(1);
  xmax=dim_lim(2);

  nbTot=size(p,2);

  cpt_boundary=0;
  cpt_interior=0;
  for i=1:nbTot
      if ((p(1,i)==xmin)|(p(1,i)==xmax)|(p(2,i)==xmin)|(p(2,i)==xmax))
	  cpt_boundary=cpt_boundary+1;
	  ind_boundary(cpt_boundary)=i;        
      else
	  cpt_interior=cpt_interior+1;
	  ind_interior(cpt_interior)=i;  
      end
  end
  
end
