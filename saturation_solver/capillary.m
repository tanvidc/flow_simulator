%Plots fw vs s curve for capillary pressure case
%Uses Brooks-Corey correlation to get pc
%Tanvi Chheda May 2016

clear all
s=0.2:0.01:0.9;
visc_w = 1;
visc_o =10;
n  = length(s);
fw = zeros(n,1);
dfwds = zeros(n,1);
dx = 0.1;

%---Input Params 
nw   = 3.5;
no   = 2;
krwe = 0.7;
kroe = 0.6;
swc  = 0.2; sor = 0.1; 
SE   = (s-swc)/(1-swc-sor); %Effective sat
Krw  = krwe*(SE).^nw;    %relative perm of w
Kro  = kroe*(1-SE).^no;

pc = 75*SE.^(-1/2);
% plot(s,pc,'Marker','*');

% dp = zeros(n); a=zeros(n);
% for i=2:n-1
%     a(i)=pc(i+1)-pc(i);
%     dp(i)=(pc(i+1)-pc(i))/dx;
% end
dp=zeros(n,1);
dp(1)=0; dp(2)=0;
dp(3:n-1)=(pc(4:n)-pc(2:n-2))./(2*dx);
dp(n)=(pc(n)-pc(n-1))/dx;

Lam_w = Krw/visc_w;
Lam_o = Kro/visc_o;
Lam_t = Lam_w+Lam_o;

dLam_w_ds = krwe*nw/(visc_w*(1-swc-sor))*((s-swc)/(1-swc-sor)).^(nw-1);
dLam_o_ds = -kroe*no/(visc_o*(1-swc-sor))*((1-s-sor)/(1-swc-sor)).^(no-1);

dp=dp.';
fw=(Lam_w./Lam_t)-(Lam_w.*Lam_o)./Lam_t.*dp;


%(a/b)' = (a'*b-b'*a)/b^2
ap = dLam_w_ds; b = Lam_w+Lam_o; bp = dLam_w_ds+dLam_o_ds; a = Lam_w;
dfwds = (ap(:).*b(:)-bp(:).*a(:))./b(:).^2;

% figure(1); hold on;
% plot(s,Krw,s,Kro,'LineWidth',2);
% legend('Krw','Kro')
% ylim([0,1]); xlabel('Sw'); 
% ylabel('Relative Permeability')

figure(2); hold on;
[fwp,~]=frac(s,1,10);
plot(s,fw,s,fwp,'LineWidth',2)
xlabel('s'); ylabel('fw');
legend('with Pc','without Pc', 'Location','southeast')

% figure(3)
% plot(s,dfwds,'LineWidth',2)
% xlabel('sw'); ylabel('dfw/ds')
