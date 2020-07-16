%Compressible Solver
%Solves the discretized form of  d(rho*phi)/dt -d/dx (lambda * dp/dx) = q
%Tanvi Chheda
%9 May 2016

clear all 

%Define parameters, interfaces
%|---X---|---X---|...|---X---|
Dirichelt_BC = 0; pL=1;pR=0;
AddWells = 1;

L = 1.0; % Length [m]
N = 100; % # grid cells
Dx = L / N; % grid size
x = linspace(Dx/2,L-Dx/2, N); %Grid centers
xi = linspace(0, L, N + 1); %interfaces

%Temporal
Dt=0.01; %time step
nTimeSteps=40;
maxNewtonIter = 50; %max Newton-loops

cf = 1; %effective compressibility
pinit=zeros(N,1); %initial pressure

%Single phase convective co-eff 
lambda = ones(N,1); % constructs vector
%lambda=1+0.3*randn(N,1); %Random heterogenous
phi = 0.4;
rho0 = 1;
p0 = 0;

%% Pressure Solver

p = zeros(N,nTimeSteps);
p(:,1)=pinit;

figure
plot(x,p)
hold on

% Repeat for each time step
for step = 2:nTimeSteps
    %1-pressure at time step n
    p_n = p(:,step-1);
    %2-initial guess for iteration loop, nu
    p_nu = p_n;
    %3- update quantities at time step n
    [rho_n,~] = density(p0, rho0, cf, p_n);

    for it = 1:maxNewtonIter
        A = zeros(N,N); %convective matrix
        C = zeros(N,N); 
        q = zeros(N,1);
        
        % 4- update properties
        [rho_nu, drhodp_nu] = density(p0, rho0, cf, p_nu);
        % 5- Calculate Residual
        T = trans(Dx,N,lambda,rho_nu);
        % 6- Construct C and A in (C+A)dp=R
        for i=1:N
            if (i>1)
                A(i,i-1) = -T(i);
                A(i,i) = T(i);
            end
            if (i<N)
                A(i,i+1) = -T(i+1);
                A(i,i) = A(i,i) + T(i+1);
            end
            C(i,i) = 1/(Dt)*phi*drhodp_nu; 
        end
        % implement BC
        
        if Dirichelt_BC==1
            i=1;
            %T(1)*p(1) + T(2)*(p(1)-p(2)) = q(1) + T(1)*pL
            A(i,i)=A(i,i)+T(1);
            q(i)=q(i)+T(1)*pL;
            %T(N)*(p(N)-p(N-1) + T(N+1)*p(N) = q(N) + T(N+1)*pR
            i=N;
            A(i,i)=A(i,i)+T(N+1); %????
            q(i)= q(i) + T(N+1)*pR;
        end
        
        if AddWells==1
            PIw1=1000; PIw2=1000; pw1=1; pw2=0;
            A(1,1)= A(1,1)+ PIw1*lambda(1);
            q(1) = q(1)+ PIw1*pw1*lambda(1);
            A(N,N)= A(N,N)+ PIw2*lambda(N);
            q(N)= q(N)+PIw2*pw2*lambda(N);
        end
        
        % 7- calculate R
        r_nu = q - phi/Dt*(rho_nu-rho_n) -A*p_nu;
        
        % 8- solve dp=(C+A)\R
        dp = (C+A)\r_nu;
        
        % 9- p_nu_plus_1 = p_nu + dp
        p_nu_plus_1 = p_nu + dp;
        
        % 10- check convergence: max(abs(dp(:))<eps
        % if converged, store soln & exit
        if max(abs(dp(:)))<1e-4 
            p(:,step) = p_nu_plus_1;
            break %terminates innermost for loop
        else
            p_nu = p_nu_plus_1;
        end
    end
    %check if you're out due to convergence, not max allowed iterations
    if it==maxNewtonIter
        print('reached max Newton It, may have not yet converged');    
    end

end

plot(x,p)
xlabel('x'); ylabel('Pressure')

%PLOT VELOCITY
V=zeros(N,nTimeSteps);
lambdaH = zeros(N+1,1); % Initialize Lambda Avg., calc. later
lambdaH(2:N) = (2*lambda(1:N-1).*lambda(2:N))./(lambda(1:N-1)+lambda(2:N));
lambdaH(1) = lambda(1); lambdaH(N+1) = lambda(N);
for t=1:nTimeSteps
    for i = 2:N
        V(i,t)=-lambdaH(i)*(p(i,t)-p(i-1,t))/Dx;
    end
    V(1,t)=-2*lambda(1)*(p(1,t)-pL)/Dx;
    V(N+1,t)=-2*lambda(N)*(pR-p(N,t))/Dx;
end


figure;
plot(xi,V)
ylim([0,5])
ylabel('Velocity');
xlabel('x');

figure;
rov=zeros(N,1);
for i=1:N
    rov(i)=rho_n(i).*((2*V(i,nTimeSteps).*V(i+1,nTimeSteps))./(V(i,nTimeSteps)+V(i+1,nTimeSteps)));
end
plot(x,rov);
ylim([1,2])
ylabel('Mass flux = {\rho} u');