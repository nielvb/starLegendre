function [A,b,c,P_c] = HBVM_Butcher(deg, r_quad)
% Computes the Butcher tableau for HBVM by Brugnano, HBVM(deg,r_quad)
%   Methods for solving ODE y'(t)=f(y(t))
%INPUT
%   deg = degree of polynomial used to approximate y(t)
%   r_quad = nb of nodes in quadrature rule
%OUTPUT
%   A = center of tableau
%   b = bottom row of tableau
%   c = first column of tableau
%NOTE:
%   Requires that chebfun is loaded
% This code is a slightly rewritten version of the 'notree.m' function
% available on the website:
% https://people.dimai.unifi.it/brugnano/LIMbook/software.html

[c,b] = legpts(deg,[0,1]);


%% Setup the parameters of the method
P_c = [];
int_c = [];
for j=0:r_quad-1
    pj = legpoly(j,[0,1],'norm');
    P_c = [P_c,pj(c(:))];
    int_c = [int_c, int_pl(j,c(:))];
end


A = int_c*P_c'*diag(b);

end

function I = int_pl(l,c)
% computes $\int_0^c p_l(t) dt$, with p_l the l-th degree normalized
% shifted Legendre polynomial

if l==0
    I = c;
else
    pnext = legpoly(l+1,[0,1],'norm');
    pprev = legpoly(l-1,[0,1],'norm');

    I = 1/(2*sqrt(2*l+1))*(1/sqrt(2*l+3)*pnext(c)-1/sqrt(2*l-1)*pprev(c));
end
end