% 1D incompressible saturation Transport Solver IMPLICIT
% Solves 1D, 2-phase, incomp saturation equation
% Tanvi Chheda, 17 May 2016

% grid config
% |--x--|--x--|--...x--|

clear all; close all;

%% Input params
tol= 1E-6; maxIter =100;
% Geometry
L  = 1;
N  = 100;
dx = L/N;
x  = linspace(dx/2,L-dx/2,N);
xi = linspace(0,L,N+1);
%Rock props;
phi = 0.2;
k   = 1;
swc = 0.2; sor = 0.1;
visc_w = 1E-3; 
visc_o = 10E-3;

%velocit distribution, should come from Pressure Solver!
Ut    = ones(N+1,1);
sinit = swc*ones(N,1); sinit(1) = 1-sor;

%TIME STEP

SimTime = 0.069;
nt      = 1000;
Dt      = SimTime/nt;
s       = zeros(N,nt); 
s(:,1)  = sinit;

figure(1); hold on;
for t = 1:nt
    if (mod(t,50)==1) 
        plot(x,s(:,t)); end
    
    nu = 1;
    S_nu = s(:,t);
    dS = 10*tol*ones(N,1);
    S_nu_plus_1=zeros(N,1);
    
    while norm(dS)>tol && nu<maxIter
        %for each Newton Raphson Iteration
        [f,dfds]=frac(S_nu, visc_w, visc_o);
        R = zeros(N,1);
        for i=2:N
            R(i)= S_nu(i)-s(i,t)+Dt/(phi*dx)*(f(i)*Ut(i+1)-f(i-1)*Ut(i));
        end
        %Compute A or Jacobian Matrix
        A = zeros(N,N); 
        A(1,1)=1;
        for i=2:N
            A(i,i-1)=-Dt/(phi*dx)*dfds(i-1)*Ut(i);
            A(i,i)=1+Dt/(phi*dx)*dfds(i)*Ut(i+1);
        end
        dS = A\(-R);
        S_nu_plus_1 = S_nu+dS;
        nu = nu+1;
        S_nu = S_nu_plus_1;
    end
    s(:,t+1)=S_nu_plus_1;
end
xlabel('x')
ylabel('Saturation')