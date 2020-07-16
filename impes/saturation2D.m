function [ Sat ] = saturation2D( Lx, Nx, Ly, Ny, K, phi, V_x, V_y, DT, Sn, Fluid)
% 1D incompressible saturation Transport Solver
% Solves 1D, 2-phase, incomp saturation equation
% Tanvi Chheda, 17 May 2016


dx = Lx/Nx;
dy = Ly/Ny;

%sor = Fluid(2); %swc = Fluid(1);

Sat=zeros(Nx,Ny);
% Sat(1) = 1-sor;

   [fn,~] = frac(Sn,K,Fluid);
   
   for i=1:Nx
       for j=1:Ny
            I = getindex(i,j, Nx,Ny);
            Iw = getindex(i-1,j, Nx,Ny);
            Is = getindex(i,j-1, Nx,Ny);
            sx=0; sy=0;
            if i>1;
                sx = (fn(I)*V_x(i+1,j)-fn(Iw)*V_x(i,j))/dx;
            end            
            if j>1;
                sy = (fn(I)*V_y(i,j+1)-fn(Is)*V_y(i,j))/dy;
            end
            Sat(I) = Sn(I)-(DT/phi)*(sx+sy);
       end
   end
