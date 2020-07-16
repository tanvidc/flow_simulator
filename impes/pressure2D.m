function [ pr, V_x, V_y ] = pressure2D( Lx, Nx, Ly, Ny, pL, pR, lambda )
%2D Incompressible P Solver
%Solves the discretized form of -d/dx (lambda * dp/dx) = q in 2D
%Tanvi Chheda
%19 Apr 2-16

%Define parameters, interfaces


%% Grid info
dx = Lx / Nx; % grid size
dy = Ly / Ny; 

x = zeros(Nx,Ny);
y = zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        x(i,j) = dx/2 + (i-1)*dx;
        y(i,j) = dy/2 + (j-1)*dy;
    end
end

xi_v = zeros(Nx+1,Ny);
yi_v = zeros(Nx+1,Ny);
xi_h = zeros(Nx,Ny+1);
yi_h = zeros(Nx,Ny+1);

for i = 1:Nx+1
    for j = 1:Ny
        xi_v(i,j) = (i-1)*dx;
        yi_v(i,j) = dy/2 + (j-1)*dy;
    end
end
for i = 1:Nx
    for j = 1:Ny+1
        xi_h(i,j) = dx/2 +(i-1) *dx;
        yi_h(i,j) = (j-1)*dy;
    end
end

%% Mobility

lambdaH_vert = ones(Nx+1,Ny);
lambdaH_vert(1,:)=lambda(1,:);
lambdaH_vert(Nx+1,:)=lambda(Nx,:);

for i=2:Nx
    for j=1:Ny
        lambdaH_vert(i,j) = 2*lambda(i,j)*lambda(i-1,j)/(lambda(i,j)+lambda(i-1,j));
    end
end

lambdaH_hori = ones(Nx,Ny+1);
lambdaH_hori(:,1)=lambda(:,1);
lambdaH_hori(:,Ny+1)=lambda(:,Ny);

for i=1:Nx
    for j=2:Ny
        lambdaH_hori(i,j) = 2*lambda(i,j)*lambda(i,j-1)/(lambda(i,j)+lambda(i,j-1));
    end
end


%% ----------------------------------------------------------------------------

T_y = zeros (Nx+1,Ny);% Transmissibility at interfaces
T_x = zeros (Nx,Ny+1);

T_y(2:Nx,:) = lambdaH_vert(2:Nx,:)/(dx^2);
T_y(1,:)=2*lambdaH_vert(1,:)/(dx^2); 
T_y(Nx+1,:)=2*lambdaH_vert(Nx+1,:)/(dx^2);

T_x(:,2:Ny) = lambdaH_hori(:,2:Ny)/(dy^2);
T_x(:,1)=2*lambdaH_hori(:,1)/(dy^2); 
T_x(:,Ny+1)=2*lambdaH_hori(:,Ny+1)/(dy^2);


%% A-matrix

q = zeros(Nx*Ny,1); 
A = zeros(Nx*Ny,Nx*Ny); 
for i=1:Nx
    for j=1:Ny
         
        ind = getindex(i,j,Nx,Ny);
        if(i > 1) % west 
            A(ind, ind-1) = -T_y(i,j);
            A(ind,ind) = A(ind,ind) + T_y(i,j);
        end
        if(i < Nx) % east 
            A(ind, ind+1 ) = -T_y(i+1,j);
            A(ind,ind) = A(ind,ind) + T_y(i+1,j);
        end
        if(j > 1) % south 
            
            A(ind, ind-Nx ) = -T_x(i,j);
            A(ind,ind) = A(ind,ind) + T_x(i,j);
        end
        if(j < Ny) % north 

            A(ind, ind+Nx ) = -T_x(i,j+1);
            A(ind,ind) = A(ind,ind) + T_x(i,j+1);
        end
    end
end

%% BC


A(1,1)= A(1,1)+ T_y(1,1);
q(1) = q(1)+ T_y(1,1)*pL;

A(Nx*Ny,Nx*Ny)= A(Nx*Ny,Nx*Ny)+ T_y(Nx+1,Ny);
q(Nx*Ny) = q(Nx*Ny)+T_y(Nx+1,Ny)*pR;

p = A\q;

pr = zeros(Nx,Ny);
for i=1:Nx*Ny
    [a,b]=getcoordinate(i,Nx,Ny);
    pr(a,b)=p(i);
end

        
%PLOT VELOCITY
V_x = zeros(Nx+1,Ny);
for i = 2:Nx
    for j = 1:Ny
        V_x(i,j)=-lambdaH_hori(i,j)*(pr(i,j)-pr(i-1,j))/dx;
    end
end

% for j = 2:Ny
%     V_x(i,j)=-lambdaH_hori(i,j)*(pr(i,j)-pr(i-1,j))/dx;
    
V_y = zeros(Nx,Ny+1);
for i = 1:Nx
    for j = 2:Ny
        V_y(i,j)=-lambdaH_vert(i,j)*(pr(i,j)-pr(i,j-1))/dy;
    end
end


end

