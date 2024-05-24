function [cc] = StarLegendre_matrix(H,sizeH, M, v)
%StarLegendre_matrix Computes approximation to matrix ODE 
%                U'(t)v = H(t)u(t)v
%   The *-Legendre numerical procedure for solving a scalar ODE
%   U'(t) = H(t)U(t) on [-1,1], with initial value U(-1) = I.
%NOTE: interval [-1,1] is used because the numerical procedure is based on
%       Legendre polynomials, we choose Legendre polynomials on [-1,1]
%       since these are usually used in publications related to the MATLAB
%       package chebfun. The code below relies on chebfun for computing
%       Legendre coefficients.
%
%REQUIRES: chebfun, see https://www.chebfun.org/
%INPUT:
%   H = function handle of given Hamiltonian H(t,k,l), where k selects the
%       kth row and l the lth column.
%       This Hamiltonian defines the ODE U'(t) = H(t)U(t) on [-1,1] with 
%       initial conditions U(-1)=I. Hence, this needs rescaling of the ODE 
%       from the usual interval [0,T] to the interval [-1,1].
%   sizeH = size of the Hamiltonian
%   M = size of the basis of Legendre polynomials used to solve the problem
%   v = starting vector of the system
%OUTPUT:
%   cc = Mxm matrix, where column k contains the Legendre coefficients of
%       the approximation to the kth function in the solution
% Literature:
% F
%   [1] S. Pozza, and N. Van Buggenhout (2023), "A new matrix equation 
%       expression for the solution of non-autonomous linear systems of 
%       ODEs". Proc. Appl. Math. Mech., 22: e202200117. 
%       URL: https://doi.org/10.1002/pamm.202200117
%   [2] C. Bonhomme, S. Pozza, and N. Van Buggenhout (2023), "A new fast 
%       numerical method for the generalized Rosen-Zener model", ArXiv.
%       URL: https://doi.org/10.48550/arXiv.2311.04144        

%% Build the system of equations that follows from the theory of *-product
% Construct right hand side for system in star-framework
for i=1:M
    b(i,1) = (-1)^(i-1)*sqrt((2*(i-1)+1)/2); % Normalization vector
end
rhs = kron(v,b);

% Identity (block) matrix
I = kron(speye(sizeH),speye(M));

% Coefficient matrix of H(t)*heaviside(t-s) in Legendre basis
for j=1:sizeH
    for k = 1:sizeH
        % Construct coefficient matrix element (k,l) of Hamiltonian
        Ajk = genCoeffMatrix(@(t)H(t,j,k),M);
        % Plug coefficient matrix into block matrix and tensor
        AM(1+M*(j-1):M*j,1+M*(k-1):M*k) = Ajk;
    end
end




% Truncation
% ----TODO-----

% Solve the linear system
xdot = (I-AM)\kron(v,b); % Solve system

%% Compute Legendre coefficients of the solution u(t)
% Block coefficient matrix for Heaviside function
H = genCoeffMatrix(@(t) ones(size(t)),M);
HM = kron(eye(sizeH),H);

% Truncation
% -----TODO-----

% Compute Legendre coefficients
cc_vec = kron(eye(sizeH),H)*xdot; 
cc = [];
for k = 1:sizeH
    cc = [cc,cc_vec((k-1)*M+1:(k*M))];
end


end