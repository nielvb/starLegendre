% Code of the numerical experiments in Section 3 and 4 of paper [1] for the
% polynomial problem. This shows the structure of the coefficient matrix and the
% behavior of the computed (approximate) Legendre coefficients of the
% solution.
% [1] https://arxiv.org/pdf/2303.11284.pdf (to be published in ETNA)
clearvars
close all


%% Problem = example from ITVOLT, transformed from [0,tend] to [-1,1]
tend = 4;
f = @(t) -1i*(t+1)*(tend/2)^2; % given
phi = @(t) exp(-1i*((t+1)*(tend/2)).^2/2); % solution
coeffsSol = cheb2leg(chebcoeffs(chebfun(phi)),'normalized'); % Legendre coeffs solution

%% Approximation
m =100;
 
for i=1:m
    b(i,1) = (-1)^(i-1)*sqrt((2*(i-1)+1)/2); % Normalization vector
end

I = eye(m);

F = genCoeffMatrix(@(t) f(t),m); % Coefficient matrix
H = genCoeffMatrix(@(t) ones(size(t)),m);% Coefficient matrix Heaviside
Hsol = genCoeffMatrix(@(t) ones(size(t)),length(coeffsSol));% Coefficient matrix Heaviside
%%
% Asign submatrices
M = 50;
A = (I-F); A = A(1:M,1:M);
B = (I-F); B = B(1:M,M+1:m);
C = (I-F); C = C(M+1:m,1:M);
D = (I-F); D = D(M+1:m,M+1:m);
u = b(1:M);
v = b(M+1:m);

x = (I-F)\b; % Solution of big problem (infinity in theory)
xM = x(1:M);
xdot = A\u; % approximation via submatrix
HM = genCoeffMatrix(@(t) ones(size(t)),M);% Coefficient matrix Heaviside


acc = 10^-15;

ind_trunc = min(find(abs(F(1:M,M))/(max(abs(F(1:M,M))))>=acc));
A_trunc = eye(M)-[F(1:ind_trunc,1:M);zeros(M-ind_trunc,M)];
xdot_trunc = A_trunc\u;

xx = (eye(M)-A\B/D*C)\(A\u)-(eye(M)-A\B/D*C)\(A\(B*(D\v))); % Correct first M entries of x

tol = 10^-14;

%% Structure, theory and practice
BD = FE(B,tol)*FE(inv(D),tol);
fig2 = figure;
% Theory
subplot(2,3,1)
ABDC = FE(inv(A),tol)*FE(BD,tol)*FE(C,tol);
term1 = [FE(inv(eye(M)-ABDC)-eye(M),tol)];
spy(FE(term1,tol))
colterm1 = (abs(term1(:,end))/max(max(abs(term1))))>tol; colterm1 = size(term1,1)-min(find(colterm1==1));
rowterm1 = (abs(term1(end,:))/max(max(abs(term1))))>tol; rowterm1 = size(term1,2)-min(find(rowterm1==1));
legend hide
title(['$(I-ABDC)^{-1}-I$ Col band =',num2str(colterm1),'--Row band=',num2str(rowterm1)],Interpreter="latex")


subplot(2,3,2)
ABD = FE(inv(A),tol)*FE(BD,tol);
term2 = FE(inv(eye(size(ABDC))-ABDC),tol)*ABD;
spy(FE(term2,tol))
colterm2 = (abs(term2(:,1))/max(max(abs(term2))))>tol; colterm2 = size(term2,1)-min(find(colterm2==1));
rowterm2 = (abs(term2(end,:))/max(max(abs(term2))))>tol; rowterm2 = max(find(rowterm2==1));
legend hide
title(['$(I-ABDC)^{-1}ABD$ Col band =',num2str(colterm2),'--Row band=',num2str(rowterm2)],Interpreter="latex")


subplot(2,3,3)
error_predicted = abs(FE(term1*xdot-term2*b(M+1:end),tol));
spy(error_predicted)
nbnonZeros = sum(error_predicted>tol);
legend hide
title(["nonzeros=",num2str(nbnonZeros)])

% Practice
subplot(2,3,4)
ABDC_prac = A\(B/D)*C;
term1_prac = (eye(M)-ABDC_prac)\eye(M)-eye(M);
c = contourf(log10(abs(term1_prac)),[-16:0]);
colorbar
set(gca, 'YDir', 'reverse' )

colterm1 = (abs(term1_prac(:,end))/max(max(abs(term1_prac))))>tol; colterm1 = size(term1_prac,1)-min(find(colterm1==1));
rowterm1 = (abs(term1_prac(end,:))/max(max(abs(term1_prac))))>tol; rowterm1 = size(term1_prac,2)-min(find(rowterm1==1));
legend hide
title(['$(I-ABDC)^{-1}-I$ Col band =',num2str(colterm1),'--Row band=',num2str(rowterm1)],Interpreter="latex")


subplot(2,3,5)
ABD_prac= A\(B/D);
term2_prac = (eye(size(ABDC_prac))-ABDC_prac)\ABD_prac;
c = contourf(log10(abs(term2_prac)),[-16:0]);
colorbar
set(gca, 'YDir', 'reverse' )
colterm2 = (abs(term2_prac(:,1))/max(max(abs(term2_prac))))>tol; colterm2 = size(term2_prac,1)-min(find(colterm2==1));
rowterm2 = (abs(term2_prac(end,:))/max(max(abs(term2_prac))))>tol; rowterm2 = max(find(rowterm2==1));
legend hide
title(['$(I-ABDC)^{-1}ABD$ Col band =',num2str(colterm2),'--Row band=',num2str(rowterm2)],Interpreter="latex")


subplot(2,3,6)
error_prac = abs(abs(xx-xdot)./abs(xx));%abs(FE(term1_prac*(A\u)-term2_prac*b(M+1:end),tol));
spy(error_prac)
nbnonZeros_prac = sum(error_prac>tol);
legend hide
title(["nonzeros=",num2str(nbnonZeros_prac)])


%% error on x
figure
semilogy(abs(xx-xdot),'ro')
hold on
semilogy(abs(xx-xdot_trunc),'g+')


%% Quality of ODE solution
% Coeffs
approxCoeffs = HM*xdot_trunc;% approxCoeffs = approxCoeffs(1:end-1);
cut = max(find(abs(approxCoeffs)>eps));
cc = HM*xdot;
%cc = cc(1:cut);
%cc = leg2cheb(cc(1:cut),'norm');

%M = M-1;
if M<length(coeffsSol) 
    approxCoeffs = [approxCoeffs(:);zeros(length(coeffsSol)-M,1)];
else
    coeffsSol = [coeffsSol(:);zeros(M-length(coeffsSol),1)];
end

figure
subplot(1,2,1)
semilogy(abs(coeffsSol-approxCoeffs),'g^')
hold on
semilogy(abs(coeffsSol-[cc;zeros(M-length(cc),1)]),'ro')

subplot(1,2,2)
semilogy(abs(coeffsSol),'b*')
hold on
semilogy(abs(approxCoeffs),'g^')
semilogy(abs(cc),'ro')


% Function evaluation
chebcc = leg2cheb(approxCoeffs(1:M-1),'norm');
approx = chebfun(chebcc,'coeffs');

xeval = linspace(-1,1,2*m);
figerr = figure;
subplot(2,2,1)
plot(xeval,real(phi(xeval)),'b')
hold on
plot(xeval, real(approx(xeval)),'r--')
legend hide
subplot(2,2,2)
plot(xeval,imag(phi(xeval)),'b')
hold on
plot(xeval,imag(approx(xeval)),'r--')
legend hide
subplot(2,2,[3,4])
err = abs(phi(xeval)-approx(xeval));
plot(xeval,err,'g+')
legend hide


%% Contour
fff = figure;
 c = contourf(log10(abs(F)));
 colorbar
 set(gca, 'YDir', 'reverse' )

%%

acc = 10^-14;
fig1 = figure;
subplot(2,2,1)
FM = F(1:M,1:M);
AA = eye(M)-FM;
logAA = log10(abs(AA)); logAA(find(isinf(logAA))) = -40; % replace -inf by -40 for dat-file
c = contourf(logAA,[-4:0]);
colorbar
set(gca, 'YDir', 'reverse' )
bandF = (abs(F(1,:))/max(max(abs(F))))>tol; bandF = max(find(bandF==1));
title(["$(I-F)$ band size=",num2str(bandF)],Interpreter="latex")


% Inverses and numerical band


Dinv = inv(D); subplot(2,2,3), 
c = contourf(log10(abs(Dinv)),[-16:0]);
colorbar
set(gca, 'YDir', 'reverse' )
bandD = (abs(Dinv(1,:))/max(max(abs(Dinv))))>tol; bandD = max(find(bandD==1));
legend hide
title(["$D^{-1}$ band size=",num2str(bandD)],Interpreter="latex")

Ainv = inv(A); subplot(2,2,2), spy(FE(Ainv,tol));
c = contourf(log10(abs(Ainv)),[-16:0]);
colorbar
set(gca, 'YDir', 'reverse' )
bandA = (abs(Ainv(1,:))/max(max(abs(Ainv))))>tol; bandA = max(find(bandA==1));
legend hide
title(["$(I_M-F_M)^{-1}$ band size=",num2str(bandA)],Interpreter="latex")




% Solve system and determine truncation error

subplot(2,2,4)
err_x = abs(xM-xdot)./abs(xM);
err_x(find(err_x==0)) = tol;
semilogy(err_x,'g+')
hold on
ind = find(err_x>tol); ind = ind(1);
semilogy([ind,ind],[eps,max(err_x)],'k--','LineWidth',1.5)
xt = [0:100:M];
indticks = find(xt<= ind); indticks = indticks(end);
xt = [xt(1:indticks),ind,xt(indticks+1:end)];
xticks(xt)
ylabel('$\vert \dot{x}-x_M\vert/\vert x_M\vert$',Interpreter='latex')
legend hide
title(['inaccurate=', num2str(M-ind)])
