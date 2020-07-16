function [ fw, dfwds ] = frac( s, K, Fluid)
fw = zeros(size(s));
dfwds = zeros(size(s));

krwe   = Fluid(7); kroe   = Fluid(8); 
swc    = Fluid(1); sor    = Fluid(2);
visc_w = Fluid(3); visc_o = Fluid(4);
nw     = Fluid(5); no     = Fluid(6);


[Lam_w, Lam_o] = lamb(s, K, Fluid);


dLam_w_ds = krwe*nw/(visc_w*(1-swc-sor))*((s-swc)/(1-swc-sor)).^(nw-1);
dLam_o_ds = -kroe*no/(visc_o*(1-swc-sor))*((1-s-sor)/(1-swc-sor)).^(no-1);

fw(:) = Lam_w(:)./(Lam_w(:)+Lam_o(:));

%(a/b)' = (a'*b-b'*a)/b^2
ap = dLam_w_ds; b = Lam_w+Lam_o; bp = dLam_w_ds+dLam_o_ds; a = Lam_w;
dfwds = (ap(:).*b(:)-bp(:).*a(:))./b(:).^2;


end