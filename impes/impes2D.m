% 2D incompressible solver
% Tanvi Chheda June 2016
clear all;

Lx  = 1;
Nx  = 20;
Ly  = 1;
Ny  = 20;
pL  = 1;
pR  = 0;
phi = 0.3;
K   = 0.1;

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

%Fluid props
swc    = 0.2;  sor    = 0.1;
nw     = 3.5;  no     = 2;
krwe   = 0.7;  kroe   = 0.6;
visc_w = 1E-3; visc_o = 10E-3;

Fluid  = [swc, sor, visc_w, visc_o, nw, no, krwe, kroe];

lambda = K/phi * ones(Nx,Ny);
Sn = swc* ones(Nx,Ny); Sn(1,1)=1-sor;
[Lw,Lo] = lamb(Sn, K, Fluid);
lambda  = Lw+Lo;

SimTime = 1; Time = 0;

%DT = 1E-6; V_x = ones(Nx+1,Ny); V_y = ones(Nx,Ny+1);
%%

figure(1); hold on; i=1;
while i<15000 % Time < SimTime %   
    [p, V_x, V_y]  = pressure2D(Lx, Nx, Ly, Ny, pL, pR, lambda);  

    DT      = Waqt(K, Fluid, phi, sor, swc, Nx, Ny, dx, dy, V_x,V_y);
    S       = saturation2D(Lx, Nx, Ly, Ny, K, phi, V_x, V_y, DT, Sn, Fluid);
    [Lw,Lo] = lamb(S, K, Fluid);
    lambda  = Lw+Lo;
    Sn = S;
    Time    = Time+DT;
    i=i+1;
end
surf(x,y,Sn) 
grid on