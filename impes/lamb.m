function [ Lam_w, Lam_o ] = lamb( s, K, Fluid)

krwe   = Fluid(7); kroe   = Fluid(8); 
swc    = Fluid(1); sor    = Fluid(2);
visc_w = Fluid(3); visc_o = Fluid(4);
nw     = Fluid(5); no     = Fluid(6);

SE   = (s-swc)/(1-swc-sor);
Krw    = (krwe*(SE)).^nw;    %relative perm of w
Kro    = (kroe*(1-SE)).^no;

%------

Lam_w = K.*Krw/visc_w;
Lam_o = K.*Kro/visc_o;

