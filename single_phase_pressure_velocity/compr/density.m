function [ rho, drhodp ] = density( p0, rho0, cf, p )

rho = rho0*(1+cf*(p-p0));
drhodp = rho0*cf;

end

