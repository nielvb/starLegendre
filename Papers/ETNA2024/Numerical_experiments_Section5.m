% % This file contains the three numerical experiments that are performed 
% in  https://arxiv.org/pdf/2303.11284.pdf
clearvars
close all


%% Set up problem and solution
choose = 1; 
if choose==1 % See Section 5.1 Toy problem
    % Parameters
    omega=5;
    beta = 1;
    M = 100;
    % Function and solution
    f = @(t) -1i*omega/beta*sin(omega*(t+1)); % Given function
    phi = @(t) exp(-1i/beta*(1-cos(omega*(t+1)))); % Solution to ODE
    f_HBVM = @(t) -2*1i*omega/beta*sin(2*omega*t); % function on [0,1]
    phi_HBVM = @(t) exp(-1i/beta*(1-cos(omega*(2*t-1)+omega))); % solution on [0,1]
elseif choose==2 % See Section 5.2 A polynomial problem
    % Parameters
    tend = 25; % Solve on interval [0,tend]
    M = 500;
    % Function and solution
    f = @(t) -1i*(t+1)*(tend/2)^2; 
    phi = @(t) exp(-1i*((t+1)*(tend/2)).^2/2);
    f_HBVM = @(t) -1i*(tend)^2*t; % on [0,1]
    phi_HBVM = @(t) exp(-1i*tend^2*t.^2/2);
elseif choose ==3 % See Section 5.3 NMR-inspired problem
    % Parameters
    nu = 5000;
    tend = 10^-2; % Solve on interval [0,tend]
    M = 1000;
    % Function and solution
    [H,Sol,~,alpha,beta] = NMR_example(nu,tend);
    i = 1; % Select one of the elements in the matrix valued function
    f = @(t) -1i*tend*pi*H(t,i,i); 
    phi = @(t) Sol(t,i,i);
    cte = -tend*2*1i*pi;
    alpha = alpha(i,i);
    beta = beta(i,i);
    f_HBVM = @(t) cte*(alpha+beta*cos(2*pi*nu*tend*t)+beta*cos(4*pi*nu*tend*t));
    phi_HBVM = @(t) exp(cte*alpha*t+cte*beta/(2*pi*nu*tend)*sin(2*pi*nu*tend*t)+...
                cte*beta/(4*pi*nu*tend)*sin(4*pi*nu*tend*t));
end



 
%% Approximation
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
semilogy(abs(coeffsSol),'b*','DisplayName','Exact')
hold on
semilogy(abs(cc_trunc),'g^','DisplayName','Approximate')
xlabel('$i$','Interpreter', 'Latex')
ylabel('$\vert c_i \vert$','Interpreter', 'Latex')
title('Legendre coefficients of solution')
legend hide

subplot(1,2,2)
semilogy(abs(coeffsSol-cc_trunc),'g^')
xlabel('$i$','Interpreter', 'Latex')
ylabel('$\textrm{err}_c$',Interpreter='latex')
title(['Error on Legendre coefficients'])
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
title('Error on real part')
subplot(2,2,4)
err2 = imag(phi(xeval)-approx_trunc(xeval));
plot(xeval,abs(err2),'g+')
title('Error on imaginary part')
legend hide

% Error on function evaluations
errf = norm(phi(xeval)-approx_trunc(xeval),"inf")/norm(phi(xeval),"inf")

%% Plot pseudospectra of (I-F)
fig_spec = figure(99);
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

%% Table providing different measure for the error and the metrics to determine invertibility of (I-F)
names = {'$Coefficient error$','$Function error$','sum alpha','numerical radius bound'};
Tab2 = array2table([max(errc),max(errf),sumalpha,numrad],'VariableNames',names);
disp(Tab2) 

%% sHBVM (see https://people.dimai.unifi.it/brugnano/LIMbook/software.html)
u0 = 1; % initial condition
h = 1; % step size
deg = M; % Degree of Legendre polynomial basis
r_quad = round(deg/2); % nb of nodes for Gauss-Legendre quadrature rule

[AA,b,c,Ps] = HBVM_Butcher(deg, r_quad);

fact = f_HBVM(h*c(:));
A = eye(deg)-h*diag(fact)*AA;
rhs = fact.*repmat(u0,deg,1);
k = A\rhs;

% Advance the solution
approx = u0 + h*sum(k(:).*b(:));


%% Compare error of sHBVM and our approach
disp('We compare the errors of the sHBVM approach and our approach:')
err_HBVM = abs(approx - phi_HBVM(h))
err_star = max(errf)
