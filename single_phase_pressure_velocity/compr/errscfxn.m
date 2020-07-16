function [ pr ] = errscfxn( Dt )

%Function takes Dt and returns pressure in final time step. 
%To be used in script errorcalc to plot error-Dt loglog plot
%Slightly compressible Solver. 
%Tanvi Chheda

%Define parameters, interfaces
%|---X---|---X---|...|---X---|

Dirichelt_BC = 1;
AddWells = 0;
L = 1.0; % Length [m]
N = 20; % # grid cells
Dx = L / N; % grid size
x = linspace(Dx/2,L-Dx/2, N); %Grid centers
xi = linspace(0, L, N + 1); %interfaces

Sim_time = 0.1;
nt = Sim_time/Dt;
c_eff = 1; %effective compressibility
poro = 0.4; %porosity

pL = 1; pR = 0;
pini = zeros(N,1); %initial pressure
p = zeros(N,nt); 
p(:,1)=pini;
PIw1=10; PIw2=10; pw1=1; pw2=0;

Lambda = zeros(N,1); % constructs vector
Lambda(1:N)=1;
LambdaH = zeros(N+1,1); % Initialize Lambda Avg., calc. later
LambdaH(2:N) = (2*Lambda(1:N-1).*Lambda(2:N))./(Lambda(1:N-1)+Lambda(2:N));
LambdaH(1) = Lambda(1); LambdaH(N+1) = Lambda(N);


T = zeros (N+1,1); % Transmissibility at interfaces
T(2:N) = LambdaH(2:N)/(Dx^2);
T(1) = 2*LambdaH(1)/(Dx^2); T(N+1)=2*LambdaH(N+1)/(Dx^2);

%% ----------------------------------------------------------------------------

C=zeros(N,N);
q = zeros(N,1); %can add source term later here
A = zeros(N,N);
for i=1:N
    C(i,i)=poro*c_eff/Dt;
end

%T(i)*(p(i)-p(i-1))+T(i+1)*(p(i)-p(i+1))=q(i)
for i=1:N
    if(i>1)
        A(i,i-1) = -T(i);
        A(i,i) = T(i);
    end

    if(i<N)
        A(i,i+1) = -T(i+1);
        A(i,i) = A(i,i)+ T(i+1); 
    end
end

% Boundary Condition
if Dirichelt_BC==1
    i=1;
    A(i,i) = A(i,i)+T(1);
    q(i) = q(i)+T(1)*pL;

    i=N;
    A(i,i) = A(i,i) + T(N+1);
    q(i) = q(i)+T(N+1)*pR;   
end

if AddWells==1
    A(1,1)= A(1,1)+ PIw1*Lambda(1);
    q(1) = q(1)+ PIw1*pw1*Lambda(1);
    A(N,N)= A(N,N)+ PIw2*Lambda(N);
    q(N)= q(N)+PIw2*pw2*Lambda(N);
end


for t = 1:nt
    %Implicit: (C+A)*p(n+1) = q + c*p(n)
    p(:,t+1) = (C+A) \ ( q(:) + C*p(:,t));
    % Explicit: P{n+1} = C/(q + C*P(n)-A*p(n))
    %p(:,t+1) = C\(q(:) + C*p(:,t) - A*p(:,t));
end
pr=p(:,nt+1);


