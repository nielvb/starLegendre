function [approx,cc] = StarLegendre_scalar(f, M, Legcoeffs)
%StarLegendre_scalar Computes approximation to scalar ODE u'(t) = f(t)u(t)
%   The *-Legendre numerical procedure for solving a scalar ODE
%   u'(t) = f(t)u(t) on [-1,1], with initial value u(-1) = 1.
%NOTE: interval [-1,1] is used because the numerical procedure is based on
%       Legendre polynomials, we choose Legendre polynomials on [-1,1]
%       since these are usually used in publications related to the MATLAB
%       package chebfun. The code below relies on chebfun for computing
%       Legendre coefficients.
%
%REQUIRES: chebfun, see https://www.chebfun.org/
%INPUT:
%   f = function handle of given scalar function f(t) defining the ODE 
%       u'(t) = f(t)u(t) on [-1,1] with initial condition u(-1)=1.
%       Hence, this needs rescaling of the ODE from the usual interval 
%       [0,T] to the interval [-1,1].
%   M = size of the basis of Legendre polynomials used to solve the problem
%   Legcoeffs = (optional) Legendre coefficients for f(t)
%OUTPUT:
%   approx = function handle for the approximation given in terms of a
%           truncated Chebyshev series with N<M terms.
%   cc = N<M Legendre coefficients of the approximation.
% Literature:
% Details on the numerical procedure and its numerical behavior, such as
% rate of convergence and truncation error, can be found in [1]. This
% reference discusses the scalar problem.
%   [1] S. Pozza, and N. Van Buggenhout (2023), "A new Legendre polynomial-based 
%      approach for non-autonomous linear ODEs" ArXiv.
%       URL: https://doi.org/10.48550/arXiv.2303.11284


%% Build the system of equations that follows from the theory of *-product
% Right-hand side, i.e., Legendre polynomials evaluated in t=-1
for i=1:M
    b(i,1) = (-1)^(i-1)*sqrt((2*(i-1)+1)/2); % Normalization vector
end

% Identity matrix of size M
I = eye(M);

% Coefficient matrix of f(t)*heaviside(t-s) in Legendre basis
if nargin >2
    F = genCoeffMatrix(@(t) f(t),M,Legcoeffs);
else
    F = genCoeffMatrix(@(t) f(t),M);
end

% Truncation of the linear system of equations, truncation is determined by 
% the requested accuracy acc.
ind_trunc = min(find(abs(F(1:M,M))/(max(abs(F(1:M,M))))>=10^-15));
A_trunc = eye(M)-[F(1:ind_trunc,1:M);zeros(M-ind_trunc,M)]; 

% Solve the linear system
xdot = A_trunc\b;

%% Compute Legendre coefficients of the solution u(t)
% Coefficient matrix for the heaviside function heaviside(t-s)
HM = genCoeffMatrix(@(t) ones(size(t)),M);

% Truncate the heaviside coefficient matrix
HM_trunc = HM; HM_trunc(M,:) = 0;

% Compute Legendre coefficients
cc = HM_trunc*xdot;

% The approximating Legendre series as a function handle
approx = chebfun(leg2cheb(cc,'norm'),'coeffs');

end