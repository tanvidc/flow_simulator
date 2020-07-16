function [s, maxdfds ] = shkfront( viscratio )
% 1D incompressible saturation Transport Solver
% Solves 1D, 2-phase, incomp saturation equation
% Tanvi Chheda, 17 May 2016

% grid config
% |--x--|--x--|--...x--|

close all;

%% Input params
% Geometry
L  = 1;
N  = 1000;
dx = L/N;
x  = linspace(dx/2,L-dx/2,N);
xi = linspace(0,L,N+1);
%Rock props;
phi = 0.3;
swc = 0.2; sor = 0.1;
visc_w = 1E-3; 
visc_o = visc_w*viscratio;

%velocity distribution, should come from Pressure Solver!
Ut      = 0.8*ones(N+1,1);
Ut(1)   = 0;
Ut(N+1) = 0;

%TIME STEP
%Find max of dfwds based on most restrictive Dt calculation.
CFL       = 0.5;
dsw = (1-sor-swc)/(N-1);
s_temp    = swc:dsw:1-sor;  %%CHANGE
[~,dfwds] = frac(s_temp, visc_w, visc_o);
max_Lam   = max(dfwds(:))*max(Ut(:))/phi;
DT        = CFL*dx/max_Lam;

% SimTime   = 10;
nt        = 5000;
s         = zeros(N,nt);

% solve for all saturations
s(:,1) = swc*ones(N,1); %IC- initial sw everywhere
s(1,1) = (1-sor); % Dirichlet BC, forced to be 1, 100% flooded

maxdfds=zeros(nt,1);

% plot(x,s(:,1));
for t = 1:nt
    [fn,dfn_ds] = frac(s(:,t),visc_w, visc_o);
    if (t>10)
        DT = calc_DT(CFL,phi,dx,Ut,dfn_ds,N);
    end
    [~, maxdfds(t)]=max(dfn_ds);
    s(1,t+1) = 1-sor;
    for i = 2:N-1
        s(i,t+1) = s(i,t) - DT/dx/phi*(fn(i)*Ut(i+1)-fn(i-1)*Ut(i));
    end  
    s(N,t+1) = swc; %BC

end


end

