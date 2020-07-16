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

SimTime   = 1;
nt        = round(SimTime/DT);
s         = zeros(Nx,Ny,nt);  

% solve for all saturations
s(:,:,1) = zeros(Nx,Ny,1); %IC- initial sw everywhere
% s(1,1,1) = (1); % Dirichlet BC, forced to be 1, 100% flooded

figure(1); hold on; a=surf(s(:,:,1));
for t = 1:nt
    [fn,dfn_ds] = fracA(s(:,:,t));
    s(1,1,t+1) = 1;
    for i=2:Nx
        for j=2:Ny 
            s(i,j,t+1) = s(i,j,t)-DT/dx/phi*(fn(i,j)*V_x(i+1,j)-fn(i-1,j)*V_x(i,j))...
                - DT/dy/phi*(fn(i,j)*V_y(i,j+1)-fn(i,j-1)*V_y(i,j)); 
        end
    end
    for j=2:Ny %i=1
        s(1,j,t+1) = s(1,j,t)-DT/dx/phi*(fn(1,j)*V_x(2,j))...
                - DT/dy/phi*(fn(1,j)*V_y(1,j+1)-fn(1,j-1)*V_y(1,j)); 
    end

    for i=2:Nx %j=1
        s(i,1,t+1) = s(i,1,t)-DT/dx/phi*(fn(i,1)*V_x(i+1,1)-fn(i-1,1)*V_x(i,1))...
                - DT/dy/phi*(fn(i,1)*V_y(i,2)); 
    end
    
%     s(Nx,:,t+1) = 0; %BC
%     s(:,Ny,t+1) = 0;
    delete(a)
    a=surf(x,y,s(:,:,t));
    pause(0.05)
end

xlabel('x[m]')
ylabel('y[m]')
zlabel('Saturation[-]')
grid on
hold off
