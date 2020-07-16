% 1D incompressible solver
% Tanvi Chheda June 2016
clear all;

L  = 1;
N  = 100;
dx = L/N;
pL = 1;
pR = 0;
phi = 0.3;
K  = 0.1;

%Fluid props
swc    = 0.2;  sor    = 0.1;
nw     = 3.5;  no     = 2;
krwe   = 0.7;  kroe   = 0.6;
visc_w = 1E-3; visc_o = 10E-3;

Fluid  = [swc, sor, visc_w, visc_o, nw, no, krwe, kroe];

lambda = K/phi * ones(N,1);
Sn = swc* ones(N,1);
[~,Ut] = pressure(L, N, pL, pR, lambda);

%TIME STEP
% Find max of dfwds based on most restrictive Dt calculation.
CFL       = 0.5;
dsw       = (1-sor-swc)/(N-1);
s_temp    = swc: dsw: 1-sor;  
[~,dfds]  = frac( s_temp, K, Fluid);
max_L     = max( dfds(:)) * max(Ut(:)) /phi;
DT        = CFL* dx/ max_L;

SimTime = 0.01; Time = 0;

figure(1); hold on;
while Time < SimTime
    [p, u]  = pressure(L, N, pL, pR, lambda);  
    S       = saturation(L, N, K, phi, u, DT, Sn, Fluid);
    [Lw,Lo] = lamb(S, K, Fluid);
    lambda  = Lw+Lo;
    Sn = S;
    Time    = Time+DT;
    
%     if Time+Dt>SimTime %So we don't pass SimTime
%         DT = SimeTime-Time;
%     end

    plot(p)
end