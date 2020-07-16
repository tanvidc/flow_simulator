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
Nx = 10; 
Ny = 10;
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

%velocit distribution, should come from Pressure Solver!
V_x = ones(Nx+1,Ny);
V_y = ones(Nx,Ny+1);
% Ut(1)   = 0;
% Ut(N+1) = 0;

%TIME STEP
% Find max of dfwds based on most restrictive Dt calculation.
CFL       = 0.5;
s_temp    = 0:0.01:1;  
[~,dfwds] = fracA(s_temp);
max_Lam   = max(dfwds(:))*max(V_x(:))/phi;
DT        = CFL*dx/max_Lam;

SimTime   = 10;
nt        = round(SimTime/DT);
s         = zeros(Nx*Ny,nt);  

% solve for all saturations
s(:,1) = zeros(Nx*Ny,1); %IC- initial sw everywhere
s(1,1) = (1); % Dirichlet BC, forced to be 1, 100% flooded


for t = 1:nt
    [fn,dfn_ds] = fracA(s(:,t));
    s(1,t+1) = 1;
    for i=2:Nx
        for j=2:Ny
            ind = getindex(i,j,Nx,Ny);          
            s(ind,t+1) = s(ind,t)-DT/dx/phi*(fn(ind)*V_x(i+1,j)-fn(ind-1)*V_x(i,j))...
                - DT/dy/phi*(fn(ind)*V_y(i,j+1)-fn(ind-Nx)*V_y(i,j)); 
        end
    end
    i=1;
    for j=2:Ny
        ind = getindex(i,j,Nx,Ny);
        s(ind,t+1) = s(ind,t)-DT/dx/phi*(fn(ind)*V_x(i+1,j))...
                - DT/dy/phi*(fn(ind)*V_y(i,j+1)-fn(ind-Nx)*V_y(i,j)); 
    end
    j=1;
    for i=2:Nx
        ind = getindex(i,j,Nx,Ny);
        s(ind,t+1) = s(ind,t)-DT/dx/phi*(fn(ind)*V_x(i+1,j)-fn(ind-1)*V_x(i,j))...
                - DT/dy/phi*(fn(ind)*V_y(i,j+1)); 
    end
    
    s(Nx*Ny,t+1) = 0; %BC
end

ss = zeros(Nx,Ny,nt);
for t=1:nt
    for i=1:Nx*Ny
        [a,b]=getcoordinate(i,Nx,Ny);
        ss(a,b,t)=s(i,t);
    end
end    

figure(1); hold on;
surf(x,y,ss(:,:,t)); 
xlabel('x')
ylabel('y')
zlabel('Saturation')
