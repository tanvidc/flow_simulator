%1D Incompressible P Solver
%Solves the discretized form of -d/dx (lambda * dp/dx) = q
%Tanvi Chheda
%19 Apr 2-16

clear all

%Define parameters, interfaces
%|---X---|---X---|...|---X---|
L = 1.0; % Length [m]
pL = 1; pR = 0;
N = 100; % # grid cells
Dx = L / N; % grid size
x = linspace(Dx/2,L-Dx/2, N); %Grid centers
xi = linspace(0, L, N + 1); %intefaces

Lambda = zeros(N,1); % constructs vector
Lambda(1:N)=1; %Homogenous
% for i=1:N
%     Lambda(i)= 2*i/N; %Heterogenous- lambda varying linearly w avg 1
% end

Lambda=2*rand(N,1); %Random heterogenous

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

%T(i)*(p(i)-p(i-1))+T(i+1)*(p(i)-p(i+1))=q(i)
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
%add T(1)*(p(1)-pL) to my linear system
%T(1)*(p(1)-pL) + T(2)*(p(1)-p(2)) = q(1)
%T(1)*p(1) + T(2)*(p(1)-p(2)) = q(1) + T(1)*pL

A(i,i)=A(i,i)+T(1);
q(i)=q(i)+T(1)*pL;

%T(N)*(p(N)-p(N-1) + T(N+1)*p(N) = q(N) + T(N+1)*pR
i=N;
A(i,i)=A(i,i)+T(N+1);
q(i)= q(i) + T(N+1)*pR;

p=A\q;

%% Plotting

figure; % Create figure

subplot(4,1,1)
plot(x,Lambda,'Marker','.','LineStyle','-','color',[0 0 0],'DisplayName','L'); hold on;
plot(xi,LambdaH,'Marker','*','LineStyle','-','color',[1 0 0],'DisplayName','{\lambda}^H');
legend show; ylabel('\lambda  and  {\lambda}^H');

subplot(4,1,2)
plot(xi,T)
ylabel('Transmissivity');

subplot(4,1,3)
plot(x,p)
ylabel('Pressure');

%PLOT VELOCITY
V=zeros(N+1,1);
for i = 2:N
    V(i)=-LambdaH(i)*(p(i)-p(i-1))/Dx;
end
V(1)=-2*Lambda(1)*(p(1)-pL)/Dx;
V(N+1)=-2*Lambda(N)*(pR-p(N))/Dx;
subplot(4,1,4)
plot(xi,V)
ylim([0,2]);
ylabel('Velocity');
xlabel('x');