%Compressible Solver
%Solves the discretized form of  d(rho*phi)/dt -d/dx (lambda * dp/dx) = q
%Tanvi Chheda
%9 May 2016

clear all 

%Define parameters, interfaces
%|---X---|---X---|...|---X---|
Dirichelt_BC = 1; 

L = 20.0; % Length [m]
Dx=0.1; % grid size
N = L/Dx; % # grid cells
x = linspace(Dx/2,L-Dx/2, N); %Grid centers
xi = linspace(0, L, N + 1); %interfaces

%Temporal
Dt=0.0000001; %time step
nTimeSteps=10000;
maxNewtonIter = 100; %max Newton-loops

pinit= 1E6*ones(N,1); %initial pressure

%Single phase convective co-eff 
K = 1E-12;
mu= 1E-17;
lambda = K/mu*ones(N,1); % constructs vector

phi = 0.2;
rho0 = 0.5;
p0 = 1E5;

%% Pressure Solver

p = 1E6*ones(N,nTimeSteps);
p(:,1)=pinit;


figure
hold on
% Repeat for each time step
for step = 2:nTimeSteps
    %1-pressure at time step n
    p_n = p(:,step-1);
    %2-initial guess for iteration loop, nu
    p_nu = p_n;
    %3- update quantities at time step n
    [rho_n,~] = density(p0, rho0, p_n);
    q = zeros(N,1);
    for it = 1:maxNewtonIter
        A = zeros(N,N); %convective matrix
        C = zeros(N,N); 
                
        % 4- update properties
        [rho_nu, drhodp_nu] = density(p0, rho0, p_nu);
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
            pW = 1E5;
                        
            A(100,101)=A(100,101)+T(101);
            q(100)=q(100)+T(101)*pW;
            
            A(101,100)=A(101,100)+T(101);
            q(101)= q(101) + T(101)*pW;
        end
        
        % 7- calculate R
        r_nu = q - phi/Dt*(rho_nu-rho_n) -A*p_nu;
        
        % 8- solve dp=(C+A)\R
        dp = (C+A)\r_nu;
        
        % 9- p_nu_plus_1 = p_nu + dp
        p_nu_plus_1 = p_nu + dp;
        
        % 10- check convergence: max(abs(dp(:))<eps
        % if converged, store soln & exit
        if max(abs(dp(:)))<10 
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
   if (mod(step,500)==1) 
        plot(x,p(:,step)); end
end


% plot(x,p)
xlabel('x'); ylabel('Pressure')

%PLOT VELOCITY
V=zeros(N+1,nTimeSteps);
lambdaH = zeros(N+1,1); % Initialize Lambda Avg., calc. later
lambdaH(2:N) = (2*lambda(1:N-1).*lambda(2:N))./(lambda(1:N-1)+lambda(2:N));
lambdaH(1) = lambda(1); lambdaH(N+1) = lambda(N);
figure;
hold on;
for t=1:nTimeSteps
    for i = 2:N
        V(i,t)=-lambdaH(i)*(p(i,t)-p(i-1,t))/Dx;
    end
    if (mod(t,500)==1) 
        plot(xi,V(:,t)); end
end

ylabel('Velocity');
xlabel('x');

% FLUX
% figure;
% rov=zeros(N,1);
% for i=1:N
%     rov(i)=rho_n(i).*((2*V(i,nTimeSteps).*V(i+1,nTimeSteps))./(V(i,nTimeSteps)+V(i+1,nTimeSteps)));
% end
% plot(x,rov);
% ylabel('Mass flux = {\rho} u');