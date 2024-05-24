% For details, see https://arxiv.org/abs/2311.04144

close all
clearvars

% Problem setup
nu = 70000;
nu0 = 20000;
tau= 1.1*10^-5;
tend = pi/nu0;



HH = @(t,nu,nu0,tau) tend*1i*pi*nu0*[0,sech(pi*t*tend/tau)*exp(-2*pi*1i*2*nu*t*tend);sech(pi*t*tend/tau)*exp(-2*pi*1i*2*nu*t*tend),0];
indexat = @(expr, k,l) expr(k,l);
H = @(t,k,l) indexat(HH(t,nu,nu0,tau),k,l);

v = zeros(2,1); v(1) =1; % Starting state

% Approximation via star-framework
M = 1500;
[cc] = StarLegendre_matrix(H,2, M, v);

figure
semilogy(abs(cc),'g.')


%% Compare to solution
cc_fct = cc(:,1);
cc_fct = cc_fct(1:600);

cc2_fct = cc(:,2);
cc2_fct = cc2_fct(1:600);


approx_fct = chebfun(leg2cheb(cc_fct,'norm'),'coeffs');
approx2_fct = chebfun(leg2cheb(cc2_fct,'norm'),'coeffs');

x_eval = linspace(-1,1,1000);
figure
subplot(2,1,1)
plot(x_eval,real(approx_fct(x_eval)),'b')
hold on
plot(x_eval,real(approx2_fct(x_eval)),'r--')
subplot(2,1,2)
plot(x_eval,imag(approx_fct(x_eval)),'b')
hold on
plot(x_eval,imag(approx2_fct(x_eval)),'r--')