% Example with a diagonal matrix valued function, one of the Eniko examples
% Since we are dealing with a diagonal problem, we know the analytical
% expression for the solution.
close all
clearvars

% Problem setup
nu = 500; % Speed of oscillations of MAS
tend = 10^-3; % Solve on interval [0,tend]
[H,U,m] = NMR_example(nu,tend); % Construct Hamiltonian and solution
H = @(t,k,l) -2*pi*1i*tend/2*H(t,k,l); 
sizeH = 16; % size of Hamiltonian 16x16

v = zeros(m,1); v(3) = 1; % Starting state

% Approximation via star-framework
M = 100;
[cc] = StarLegendre_matrix(H,sizeH, M, v);

figure
semilogy(abs(cc),'g.')


%% Compare to solution
cc_fct = cc(:,3);
cc_fct = cc_fct(1:65);

approx_fct = chebfun(leg2cheb(cc_fct,'norm'),'coeffs');

x_eval = linspace(-1,1,1000);
sol_fct = @(t) U(t,3,3);
figure
subplot(2,2,1)
plot(x_eval,real(sol_fct(x_eval)),'b')
hold on
plot(x_eval,real(approx_fct(x_eval)),'r--')
subplot(2,2,2)
plot(x_eval,imag(sol_fct(x_eval)),'b')
hold on
plot(x_eval,imag(approx_fct(x_eval)),'r--')
subplot(2,2,[3,4])
semilogy(abs(approx_fct(x_eval)-sol_fct(x_eval)),'g.')