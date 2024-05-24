function [H,U,m] = NMR_example(MASfrequency,texp)
% Generates the solution to example A (by Eniko) on interval [-1,1]
%INPUT:
%       MASfrequency = parameter that controls oscillations in
%                       Magic Angle Spinning-experiment
%       texp = time of experiment, [0;texp] will be mapped onto [-1,1]
%OUTPUT:
%       H = Hamiltonian appearing in ODE:  -i*2*pi*H(t) on [0,texp]
%       U = function in t and s, U(t,-1) representing analytic solution to
%       the ode U'(t,-1) = -2 pi sqrt(-1) A(t) U(t,-1)
%       m = size of H and U, i.e., mxm

m = 2^4;

Iz = [0.5 0; 0 -0.5 ];
%% Generate all Basis sets, will be used to construct solution:
s1z = genOperatorSingleSpin(4, 1, Iz);
s2z = genOperatorSingleSpin(4, 2, Iz);
s3z = genOperatorSingleSpin(4, 3, Iz);
s4z = genOperatorSingleSpin(4, 4, Iz);

s12z = genOperatorDoubleSpin(4, 1, 2, Iz);
s13z = genOperatorDoubleSpin(4, 1, 3, Iz);
s14z = genOperatorDoubleSpin(4, 1, 4, Iz);
s23z = genOperatorDoubleSpin(4, 2, 3, Iz);
s24z = genOperatorDoubleSpin(4, 2, 4, Iz);
s34z = genOperatorDoubleSpin(4, 3, 4, Iz);



%calculate the time ordered exponential: propagator = exp(-i*dt*Hamiltonian)
off1 = 0.2;
off2 = 0.1;
off3 = 0;
off4 = -0.2;




weakcc12 = 1000; %in Hz ca H-N coupling constant
weakcc13 = 2000; %in Hz ca H-C coupling constant
weakcc14 = 500; %in Hz ca H-P
weakcc23 = 1200; %in Hz ca N-C
weakcc24 = 1400; %in Hz ca N-P
weakcc34 = 800; %in Hz C-P



%% Generate Hamiltonian H and solution U
% Example A: Diagonal Matrix. 4 spins, with chemical shifts (offsets), all weakly coupled. different time dependent coupling constants %
alpha = s1z*off1 + s2z*off2 + s3z*off3 + s4z*off4;
beta = weakcc12*s12z + weakcc13*s13z+ weakcc14*s14z + weakcc23*s23z + weakcc24*s24z + weakcc34*s34z;

H = @(tau,k,l) (alpha(k,l) + beta(k,l)*cos(2*pi*MASfrequency*texp/2*(tau+1)) + beta(k,l)*cos(4*pi*MASfrequency*texp/2*(tau+1)));

argument = @(tau,k,l) alpha(k,l)*(tau+1) + ...
    beta(k,l)/(2*pi*MASfrequency*texp/2)*sin(2*pi*MASfrequency*texp/2*(tau+1)) +...
    beta(k,l)/(4*pi*MASfrequency*texp/2)*sin(4*pi*MASfrequency*texp/2*(tau+1));
U = @(tau,k,l) exp(-2*pi*1i*texp/2*argument(tau,k,l))*(k==l);


end


function ssys = genOperatorDoubleSpin(numberofspins, spinID1, spinID2, CBO)
    ssys = [0.5] ;
    for n = 1:numberofspins
        if n == spinID1 || n == spinID2
            ssys =  kron(ssys, CBO * 2);
        else
            ssys = kron(ssys, eye(2));
        end
    end
end

function ssys = genOperatorSingleSpin(numberofspins, spinID, CBO)
    ssys = [0.5];
    for n = 1:numberofspins
        if n == spinID 
            ssys =  kron(ssys, CBO * 2);
        else
            ssys = kron(ssys, eye(2));
        end
    end
end

