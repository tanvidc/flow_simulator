% 2D incompressible saturation Transport Solver
% Solves 2D, 2-phase, incomp saturation equation
% Tanvi Chheda, 24 May 2016

% grid config
% |--x--|--x--|--...x--|

clear all; close all;

%% Input params
% Geometry
Lx = 1; % Length [m]
Ly = 1;
Nx = 50; 
Ny = 50;
x = zeros(Nx,Ny);
y = zeros(Nx,Ny);
dx = Lx / Nx; % grid size
dy = Ly / Ny; 
for i=1:Nx
    for j=1:Ny
        x(i,j) = dx/2 + (i-1)*dx;
        y(i,j) = dy/2 + (j-1)*dy;
    end
end

%%
%Rock props;
phi = 0.3;
visc_w = 5E-3; 
visc_o = 1E-3;
swc = 0.2; sor = 0.1;

%velocit distribution, try to bring from Pressure Solver!
ut = ones(Nx+1,Ny+1);


%TIME STEP
% Find max of dfwds based on most restrictive Dt calculation.
CFL       = 0.5;
s_temp    = 0:0.01:1;  
[~,dfwds] = frac(s_temp, visc_w, visc_o);
max_Lam   = max(dfwds(:))*max(ut(:))/phi;
DT        = CFL*dx/max_Lam;

SimTime   = 1;
nt        = 150;
s= zeros(Nx*Ny,nt);

s(:,1) = swc; %IC- initial sw everywhere
s(1,:) = 1-sor;  % Dirichlet BC 

figure(1); hold on; a = surf(x,y,reshape(s(:,1),Nx,Ny));
for t = 1:nt 
   [fn,dfn_ds] = frac(s(:,t),visc_w,visc_o);
    
   for i=1:Nx
       for j=1:Ny
            I = getindex(i,j, Nx,Ny);
            Iw = getindex(i-1,j, Nx,Ny);
            Is = getindex(i,j-1, Nx,Ny);
            sx=0; sy=0;
            if i>1;
                sx = (fn(I)*ut(i+1)-fn(Iw)*ut(i))/dx;
            end            
            if j>1;
                sy = (fn(I)*ut(j+1)-fn(Is)*ut(j))/dy;
            end
            s(I,t+1) = s(I,t)-(DT/phi)*(sx+sy);
       end
   end
    delete(a);     
    a = surf (x,y,reshape(s(:,t+1),Nx,Ny));
    pause(0.1);   
end
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('saturation [-]');
    
