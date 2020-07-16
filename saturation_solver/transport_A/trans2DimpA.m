% 1D incompressible saturation Transport Solver IMPLICIT
% Solves 1D, 2-phase, incomp saturation equation
% Tanvi Chheda, 17 May 2016

% grid config
% |--x--|--x--|--...x--|

clear all; close all;

%% Input params
tol= 1E-6; maxIter =100;
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

%Rock props;
phi = 0.2;
k   = 1;
mu_w= 1E-3;
mu_o= 1E-2;

%velocit distribution, should come from Pressure Solver!
V_x = ones(Nx+1,Ny);
V_y = ones(Nx,Ny+1);
sinit = zeros(Nx*Ny,1); sinit(1,1) = 1;

%TIME STEP
SimTime = 1;
nt      = 1000;
Dt      = 0.001;
s       = zeros(Nx*Ny,nt); 
s(:,1)  = sinit;

for t = 1:nt
    
    nu = 1;
    S_nu = s(:,t);
    dS = 10*tol*ones(Nx*Ny,1);
    S_nu_plus_1=zeros(Nx*Ny,1);
    
    while norm(dS)>tol && nu<maxIter
        %for each Newton Raphson Iteration
        [f,dfds]=fracA(S_nu);
        R = zeros(Nx*Ny,1);
        for i=Nx+1:Nx*Ny
            [a,b]=getcoordinate(i,Nx,Ny);
            R(i)= S_nu(i)-s(i,t)+Dt/(phi*dx)*(f(i)*V_x(a+1,b)-f(i-1)*V_x(a,b))...
                +Dt/(phi*dy)*(f(i)*V_y(a,b+1)-f(i-Nx)*V_y(a,b));
        end
        
        %Compute A or Jacobian Matrix
        A = zeros(Nx*Ny,Nx*Ny); 
        A(1,1)=1;
        for i=2:Nx*Ny
            [a,b]=getcoordinate(i,Nx,Ny);
            A(i,i-1)=-Dt/(phi*dx)*dfds(i-1)*V_x(a,b);
            A(i,i)=1+Dt/(phi*dx)*dfds(i)*V_x(a+1,b)+...
                Dt/(phi*dy)*dfds(i)*V_y(a,b+1);
            if i>Nx
                A(i,i-Nx)=-Dt/(phi*dy)*dfds(i-Nx)*V_y(a,b); end
        end

        dS = A\(-R);
        S_nu_plus_1 = S_nu+dS;
        nu = nu+1;
        S_nu = S_nu_plus_1;
    end
    s(:,t+1)=S_nu_plus_1;
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