function [ fw, dfwds ] = frac( s, visc_w, visc_o )
fw = zeros(size(s));
dfwds = zeros(size(s));

%---Input Params
nw   = 3.5;
no   = 2;
krwe = 0.7;
kroe = 0.6;
swc = 0.2; sor = 0.1;

SE   = (s-swc)/(1-swc-sor);
Krw    = krwe*(SE).^nw;    %relative perm of w
Kro    = kroe*(1-SE).^no;

% figure(1); hold on;
% plot(s,Krw,s,Kro,'LineWidth',2);
% legend('Krw','Kro')
% ylim([0,1]); xlabel('Sw'); 
% ylabel('Relative Permeability')

%------

Lam_w = Krw/visc_w;
Lam_o = Kro/visc_o;

dLam_w_ds = krwe*nw/(visc_w*(1-swc-sor))*((s-swc)/(1-swc-sor)).^(nw-1);
dLam_o_ds = -kroe*no/(visc_o*(1-swc-sor))*((1-s-sor)/(1-swc-sor)).^(no-1);

fw(:) = Lam_w(:)./(Lam_w(:)+Lam_o(:));

%(a/b)' = (a'*b-b'*a)/b^2
ap = dLam_w_ds; b = Lam_w+Lam_o; bp = dLam_w_ds+dLam_o_ds; a = Lam_w;
dfwds = (ap(:).*b(:)-bp(:).*a(:))./b(:).^2;

% figure(2)
% plot(s,fw,'LineWidth',2)
% xlabel('s'); ylabel('fw')

% figure(3)
% plot(s,dfwds,'LineWidth',2)
% xlabel('sw'); ylabel('dfw/ds')


end