function [ pres, velo ] = pressure( L, N, pL, pR, Lambda )
%1D Incompressible P Solver
%Solves the discretized form of -d/dx (lambda * dp/dx) = q
%Tanvi Chheda


Dx = L / N; % grid size

LambdaH = zeros(N+1,1); % Initialize Lambda Avg., calc. later
LambdaH(2:N) = (2*Lambda(1:N-1).*Lambda(2:N))./(Lambda(1:N-1)+Lambda(2:N));
LambdaH(1) = Lambda(1); LambdaH(N+1) = Lambda(N);


T = zeros (N+1,1);% Transmissibility at interfaces
T(2:N) = LambdaH(2:N)/(Dx^2);
T(1)=2*Lambda(1)/(Dx^2); T(N+1)=2*Lambda(N)/(Dx^2);


%% ----------------------------------------------------------------------------

%Construct and fill our linear system
q = zeros(N,1); A = zeros(N,N);
for i=1:N
    q(i)=0; %can add source term later
end


for i=1:N
    %step1: T(i)*(p(i)-p(i-1))
    if i>1
        A(i,i-1) = -T(i);
        A(i,i) = T(i);
    end
    %step2: T(i+1)*(p(i)-p(i+1))
    if i<N
        A(i,i+1)=-T(i+1);
        A(i,i) = A(i,i)+T(i+1);
    end
end

%employ boundary condition
i=1;

A(i,i)=A(i,i)+T(1);
q(i)=q(i)+T(1)*pL;

%T(N)*(p(N)-p(N-1) + T(N+1)*p(N) = q(N) + T(N+1)*pR
i=N;
A(i,i)=A(i,i)+T(N+1);
q(i)= q(i) + T(N+1)*pR;

p=A\q;

pres = p;

%% VELOCITY
V=zeros(N+1,1);
for i = 2:N
    V(i)=-LambdaH(i)*(p(i)-p(i-1))/Dx;
end
V(1)=-2*Lambda(1)*(p(1)-pL)/Dx;
V(N+1)=-2*Lambda(N)*(pR-p(N))/Dx;

velo = V;

end