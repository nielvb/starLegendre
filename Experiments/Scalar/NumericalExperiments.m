% An experiment from the paper https://arxiv.org/pdf/2303.11284.pdf
% See Section 5.1 Toy problem
clearvars
close all


%% Set up problem and solution
choose = 3; 
if choose==1 % See Section 5.1 Toy problem
    % Parameters
    nu=50;
    cte = 200;
    % Function and solution
    f = @(t) -1i*nu/cte*sin(nu*(t+1)); % Given function
    phi = @(t) exp(-1i/cte*(1-cos(nu*(t+1)))); % Solution to ODE

elseif choose==2 % See Section 5.2 A polynomial problem
    tend = 25; % Solve on interval [0,tend]

    f = @(t) -1i*(t+1)*(tend/2)^2; 
    phi = @(t) exp(-1i*((t+1)*(tend/2)).^2/2);
    
elseif choose ==3 % See Section 5.3 NMR-inspired problem
    % Parameters
    nu = 5000;
    tend = 10^-2; % Solve on interval [0,tend]

    [H,Sol] = NMR_example(nu,tend);
    i = 1; % Select one of the elements in the matrix valued function
    f = @(t) -1i*tend*pi*H(t,i,i); 
    phi = @(t) Sol(t,i,i);
end


cf = chebfun(f);
abs(cf.coeffs(end))/abs(cf.coeffs(1))
sum(abs(cheb2leg(cf.coeffs,'normal')))

 
%% Approximation
M =1100; % Size of Legendre basis 
% Construct right hand side (Legendre polynomials evaluated in -1)
for i=1:M
    b(i,1) = (-1)^(i-1)*sqrt((2*(i-1)+1)/2); % Normalization vector
end
I = eye(M); % Identity matrix
F = genCoeffMatrix(@(t) f(t),M); % Coefficient matrix for f(t)*heaviside(t-s)
HM = genCoeffMatrix(@(t) ones(size(t)),M);% Coefficient matrix for Heaviside(t-s)

% Perform truncation at requested accuracy
acc = 10^-15; % Accuracy for truncation
ind_trunc = min(find(abs(F(1:M,M))/(max(abs(F(1:M,M))))>=acc)); % index of truncation

% Compute approximate Legendre coefficients via truncated linear system
% solve
A_trunc = eye(M)-[F(1:ind_trunc,1:M);zeros(M-ind_trunc,M)]; % Truncated matrix defining linear system of equations
xdot_trunc = A_trunc\b; % Solve truncated system
HM_trunc = HM; HM_trunc(M,:) = 0; % Truncated coefficient matrix for Heaviside function
cc_trunc = HM_trunc*xdot_trunc; % Apply Heaviside function (star-product) to obtain approximate Legendre coefficients for solution
approx_trunc = chebfun(leg2cheb(cc_trunc,'norm'),'coeffs'); % Approximate solution as a function

%% Quality of ODE solution - erros
% Legendre coefficients of solution (used for error computation)
coeffsSol = cheb2leg(chebcoeffs(chebfun(phi,'trunc',M)),'normalized'); 

% Plot exact and computed coefficients
figcoeffs = figure;
subplot(1,2,1)
semilogy(abs(coeffsSol),'b*')
hold on
semilogy(abs(cc_trunc),'g^')
xlabel('$i$','Interpreter', 'Latex')
ylabel('$\vert c_i \vert$','Interpreter', 'Latex')
legend hide

subplot(1,2,2)
semilogy(abs(coeffsSol-cc_trunc),'g^')
hold on
xlabel('$i$','Interpreter', 'Latex')
ylabel('$\textrm{err}_c$',Interpreter='latex')
legend hide

% Error on coefficients
errc = abs(coeffsSol-cc_trunc)/norm(coeffsSol,"inf");
maxerrc = max(errc)


% Plot exact and approximate solution, real and imaginary part
xeval = linspace(-1,1,10*M);
figerr = figure;
subplot(2,2,1)
plot(xeval,real(phi(xeval)),'b')
hold on
plot(xeval, real(approx_trunc(xeval)),'g--')
legend hide
title('Real part')
subplot(2,2,2)
plot(xeval,imag(phi(xeval)),'b')
hold on
plot(xeval, imag(approx_trunc(xeval)),'g--')
legend hide
title('Imaginary part')
subplot(2,2,3)
err = real(phi(xeval)-approx_trunc(xeval));
semilogy(xeval,abs(err),'ro')
subplot(2,2,4)
err2 = imag(phi(xeval)-approx_trunc(xeval));
plot(xeval,abs(err2),'g+')
legend hide

% Error on function evaluations
errf = norm(phi(xeval)-approx_trunc(xeval),"inf")/norm(phi(xeval),"inf")

% Does the last Legendre coefficients provide
err_est = abs(cc_trunc(ind_trunc-5:ind_trunc+1))/(abs(cc_trunc(1)))
%% Plot pseudospectra of (I-F)
fig_spec = figure(99)
try
    eigtool(eye(size(F))-F)
catch
    %do nothing
    close(fig_spec)
end


%% Studying the invertibility of (I-F) through bounds on the numerical radius

% Bound based on sum of Legendre coefficients
cf = chebfun(f);
sumalpha = sum(abs(cheb2leg(cf.coeffs,'normal')));

% Approximation of numerical radius
numrad = 0;
for theta = 0:0.01:2*pi
    numrad = max(numrad, max(abs(eigs(exp(1i*theta)*(F+F')/2,1,'largestabs'))));
end

%% Table
names = {'$Coefficient error$','$Function error$','sum alpha','numerical radius bound'};
Tab2 = array2table([max(errc),max(errf),sumalpha,numrad],'VariableNames',names);
disp(Tab2)  
