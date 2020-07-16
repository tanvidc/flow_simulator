function [ Sat ] = saturation( L, N, K, phi, Ut, DT, Sn, Fluid)
% 1D incompressible saturation Transport Solver
% Solves 1D, 2-phase, incomp saturation equation
% Tanvi Chheda, 17 May 2016


dx = L/N;
swc = Fluid(1); sor = Fluid(2);

[fn,~] = frac(Sn, K,Fluid);
s=zeros(size(Sn));

s(1) = 1-sor;

for i = 2:N-1
    s(i) = Sn(i) - DT/dx/phi*(fn(i)*Ut(i+1)-fn(i-1)*Ut(i));
end  
s(N) = swc; 
    
Sat = s;
end

