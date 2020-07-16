function [ DT ] = calc_DT(CFL,phi,dx,ut,dfn_ds,N)

    LAM=ones(N+1,1);
    % loop in all interfaces, find max Lam (Characteristic speed)
    LAM(1) = ut(1)/phi*(dfn_ds(1));
    LAM(N+1) = ut(N+1)/phi*(dfn_ds(N));
    for i = 2:N
        LAM(i) = ut(i)/phi*(dfn_ds(i)+dfn_ds(i-1))/2;
    end
    
    max_Lam = max(LAM(:));
    DT = CFL*dx/max_Lam;
    
%     if (max_Lam < 1E-6)
%         DT = DT_max;
%     else
%         DT = max(CFL*dx/max_Lam,DT_max);
%     end
    
end

