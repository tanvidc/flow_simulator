function [ T ] = trans( DX, N, lambda, rho )

T(1) = rho (1) *lambda(1)/(DX^2/2);
T(N+1) = rho (N) *lambda(N)/(DX^2/2);

for i = 2:N
    T(i) = rho(i)*lambda(i)*rho(i-1)*lambda(i-1)/(rho(i)*lambda(i)+rho(i-1)*lambda(i-1))/DX^2;
end

