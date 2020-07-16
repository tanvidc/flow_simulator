function [ DT ] = Waqt( K, Fluid, phi, sor, swc, Nx, Ny, dx, dy, V_x,V_y )

%TIME STEP
% Find max of dfwds based on most restrictive Dt calculation.
CFL       = 0.4;
dsw       = (1-sor-swc)/(Nx-1);
s_temp    = swc: dsw: 1-sor;  
[~,dfds]  = frac( s_temp, K, Fluid);

w         = max(V_x(:));  q = max(V_y(:));
max_V     = max(w,q);
max_L     = max( dfds(:)) * max_V /phi;
DT        = 0.1*CFL* dx/ max_L;

end

