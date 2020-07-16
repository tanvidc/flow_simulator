% 1D incompressible saturation Transport Solver
% Solves 1D, 2-phase, incomp saturation equation
% Tanvi Chheda, 17 May 2016

% grid config
% |--x--|--x--|--...x--|

clear all; close all;

%% Input params
% Geometry
L  = 100;
N  = 100;
dx = L/N;
x  = linspace(dx/2,L-dx/2,N);
xi = linspace(0,L,N+1);
%Rock props;
phi = 0.3;
swc = 0.2; sor = 0.1;
visc_w = 1E-3; 
visc_o = 10E-3;

%velocity distribution, should come from Pressure Solver!
Ut      = 0.8*ones(N+1,1);
Ut(1)   = 0;
Ut(N+1) = 0;

%TIME STEP
% Find max of dfwds based on most restrictive Dt calculation.
CFL       = 0.5;
dsw = (1-sor-swc)/(N-1);
s_temp    = swc:dsw:1-sor;  %%CHANGE
[~,dfwds] = frac(s_temp, visc_w, visc_o);
max_Lam   = max(dfwds(:))*max(Ut(:))/phi;
DT        = CFL*dx/max_Lam;

SimTime   = 10;
nt        = round(SimTime/DT);
s         = zeros(N,nt);

% solve for all saturations
s(:,1) = swc*ones(N,1); %IC- initial sw everywhere
s(1,1) = (1-sor); % Dirichlet BC, forced to be 1, 100% flooded

figure(1); hold on;
% plot(x,s(:,1));
for t = 1:nt
    [fn,dfn_ds] = frac(s(:,t),visc_w, visc_o);
    if (t>10)
        DT = calc_DT(CFL,phi,dx,Ut,dfn_ds,N);
    end
    
    s(1,t+1) = 1-sor;
    for i = 2:N-1
        s(i,t+1) = s(i,t) - DT/dx/phi*(fn(i)*Ut(i+1)-fn(i-1)*Ut(i));
    end  
    s(N,t+1) = swc; %BC
    if (mod(t,10)==1) 
        plot(x,s(:,t+1)); end
end
xlabel('x')
ylabel('Saturation')
% plot(x,s(:,t),'LineWidth',2)

%%

% ----------------------------------------------------
% Hadi Hajibeygi – All Rights Reserved – 29 April 2015
% ANALITYCAL SOLUTION OF BUCKLEY LEVERETT EQUATION
% ----------------------------------------------------

% parameters

visc_w = 1E-3; 
visc_o = 10E-3;
nw   = 3.5;
no   = 2;
ut   = 1;
krwe = 0.7;
kroe = 0.6;
porosity= 0.3;
% saturation
nsw = 101; % number of gridpoints
sw = swc:dsw:(1-sor); % saturation range
% relperms + mobility + fracflow
krw = krwe*((sw-swc)/(1-swc-sor)).^nw;
kro = kroe*((1-sw-sor)/(1-swc-sor)).^no;
mob_w = krw/visc_w;
mob_o = kro/visc_o;
fw = mob_w./(mob_w+mob_o);
% fracflow derivative (numerical)
rr = 2:(N-1);
dfw1 = [0 (fw(rr+1)-fw(rr-1))./(sw(rr+1)-sw(rr-1)) 0];
vw = ut/porosity*dfw1;
% jump velocity
dfw2 = (fw - fw(1))./(sw-sw(1));
[shock_vel, shock_index] = max(dfw2);
shock_sat = sw(shock_index);
rr2 = shock_index:N;
% valid saturation range within BL solution
sw2 = [sw(1) sw(1) sw(rr2)];
vw2 = vw(rr2);
t = 10;
vw2 = [L/t vw2(1)+10E-4 vw2];

x2 = vw2*t;
% plot of saturation as function of position
%plot(x2,sw2,'.-','linewidth',2,'color','k');

hold off
vq=interp1(x2,sw2,x);


