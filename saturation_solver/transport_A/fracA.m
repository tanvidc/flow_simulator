function [ fw, dfwds ] = frac( s )
fw = zeros(size(s));
dfwds = zeros(size(s));

%---Input Params

Krw    = s.^2;    %relative perm of w
Kro    = (1-s).^2;
visc_w = 10E-3; 
visc_o = 1E-3;

%------

Lam_w = Krw/visc_w;
Lam_o = Kro/visc_o;

dLam_w_ds = 2*s/visc_w;
dLam_o_ds = -2*(1-s)/visc_o;

fw(:) = Lam_w(:)./(Lam_w(:)+Lam_o(:));

%(a/b)' = (a'*b-b'*a)/b^2
ap = dLam_w_ds; b = Lam_w+Lam_o; bp = dLam_w_ds+dLam_o_ds; a = Lam_w;
dfwds = (ap(:).*b(:)-bp(:).*a(:))./b(:).^2;

end


