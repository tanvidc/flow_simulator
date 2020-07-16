%2D Incompressible P Solver
%Solves the discretized form of -d/dx (lambda * dp/dx) = q in 2D
%Tanvi Chheda
%19 Apr 2-16

clear all
%Define parameters, interfaces
Lx = 1; % Length [m]
Ly = 1;
Nx = 10; 
Ny = 10;
lambda = ones(Nx,Ny); 
% for i=1:Nx          % for heterogenous reservoir
%     for j= 1:Ny     % say permeability increases away from origin
%         lambda(Nx,Ny)=i+j; 
%     end
% end; 

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

% figure;
% surf(xi_v,yi_v,T_y);
% xlabel('x'); ylabel('y'),zlabel('Vertical Transmissivity')

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

%%
Dirichelt_BC = 1; pL=1; pR=0;
AddWells = 0;
if Dirichelt_BC == 1;
   co_p1 = [1,1]; 
   I_p1 = getindex(co_p1(1),co_p1(2),Nx,Ny);
   %T(1)*p(1) + T(2)*(p(1)-p(2)) = q(1) + T(1)*pL
   A(I_p1,I_p1)=A(I_p1,I_p1)+T_x(I_p1)+T_y(I_p1);
   q(I_p1)=q(I_p1)+T_x(I_p1)*pL+T_y(I_p1)*pL;
   %T(N)*(p(N)-p(N-1) + T(N+1)*p(N) = q(N) + T(N+1)*pR
   co_p2 = [Nx,Ny];
   I_p2 = getindex(co_p2(1),co_p2(2),Nx,Ny);
   A(I_p2,I_p2)=A(I_p2,I_p2)+T_x(I_p2)+T_y(I_p2);
   q(I_p2)=q(I_p2)+T_x(I_p2)*pR+T_y(I_p2)*pR;
end
%% Insert 2 wells
if AddWells == 1;
    PI_w1=1000; pw_1=1; perf_w1 = [1,1]; 
    I_w1 = getindex(perf_w1(1),perf_w1(2),Nx,Ny);

    A(I_w1,I_w1)= A(I_w1,I_w1)+ PI_w1*lambda(perf_w1(1),perf_w1(2));
    q(I_w1) = q(I_w1)+ PI_w1 * pw_1 * lambda(perf_w1(1),perf_w1(2));

    PI_w2=1000; pw_2=0; perf_w2 = [Nx,Ny]; 
    I_w2 = getindex(perf_w2(1),perf_w2(2),Nx,Ny);

    A(I_w2,I_w2)= A(I_w2,I_w2)+ PI_w2*lambda(perf_w2(1),perf_w2(2));
    q(I_w2) = q(I_w2)+ PI_w2 * pw_2 * lambda(perf_w2(1),perf_w2(2));
end

p = A\q;

pr = zeros(Nx,Ny);
for i=1:Nx*Ny
    [a,b]=getcoordinate(i,Nx,Ny);
    pr(a,b)=p(i);
end

figure;
surf(x,y,pr)
xlabel('x'); ylabel('y'),zlabel('Pressure')
        
%PLOT VELOCITY
V_x = zeros(Nx+1,Ny);
for i = 2:Nx
    for j = 1:Ny
        V_x(i,j)=lambdaH_hori(i,j)*(pr(i,j)-pr(i-1,j))/dx;
    end
end
V_y = zeros(Nx,Ny+1);
for i = 1:Nx
    for j = 2:Ny
        V_y(i,j)=lambdaH_vert(i,j)*(pr(i,j)-pr(i,j-1))/dy;
    end
end


figure
surf(xi_v,yi_v,V_x)
ylabel('y'); xlabel('xi'); zlabel('y Velocity');

figure
surf(xi_h,yi_h,V_y)
ylabel('yi'); xlabel('x'); zlabel('x Velocity');