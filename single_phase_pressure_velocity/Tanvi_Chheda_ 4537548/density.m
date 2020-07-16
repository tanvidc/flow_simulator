function [ rho, drhodp ] = density( p0, rho0, p )

rho = rho0*p/p0/0.8;
drhodp = p0*0.8/rho0;

end

