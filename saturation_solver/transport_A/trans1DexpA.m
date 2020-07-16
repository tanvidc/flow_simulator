% 1D incompressible saturation Transport Solver
% Solves 1D, 2-phase, incomp saturation equation
% Tanvi Chheda, 17 May 2016

% grid config
% |--x--|--x--|--...x--|

clear all; close all;

%% Input params
% Geometry
L  = 1;
N  = 100;
dx = L/N;
x  = linspace(dx/2,L-dx/2,N);
xi = linspace(0,L,N+1);
%Rock props;
phi = 0.2;

%velocit distribution, should come from Pressure Solver!
ut      = ones(N+1,1);
ut(1)   = 0;
ut(N+1) = 0;

%TIME STEP
% Find max of dfwds based on most restrictive Dt calculation.
CFL       = 0.5;
s_temp    = 0:0.01:1;
[~,dfwds] = frac(s_temp);
max_Lam   = max(dfwds(:))*max(ut(:))/phi;
DT        = CFL*dx/max_Lam;

SimTime   = 1;
nt        = 200; %round(SimTime/DT);
s         = zeros(N,nt);

% solve for all saturations
s(:,1) = zeros(N,1); %IC- initial sw everywhere
s(1,1) = 1; % Dirichlet BC, forced to be 1, 100% flooded

figure(1); hold on;
plot(x,s(:,1));
for t = 1:nt
    [fn,dfn_ds] = frac(s(:,t));
    if (t>10)
        DT = calc_DT(CFL,phi,dx,ut,dfn_ds,N);
    end
    s(1,t+1) = 1;
    for i = 2:N-1
        s(i,t+1) = s(i,t) - DT/dx/phi*(fn(i)*ut(i+1)-fn(i-1)*ut(i));
    end  
    s(N,t+1) = 0; %BC
    if (mod(t,10)==1) 
        plot(x,s(:,t+1)); end
end
xlabel('x')
ylabel('Saturation')


