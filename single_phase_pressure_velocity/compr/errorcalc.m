%Plots error between two sub-sequent time step refinements
%solutions obtained by halving the time step m-times
%log-log plot, slope gives the order of consistency!

m=7; Dt=zeros(1,m); Dt(1)=0.1; e=zeros(1,m); x=0:m-1; 
for i=2:m
    Dt(i)=Dt(i-1)/2;
    A=errscfxn(Dt(i-1));
    B=errscfxn(Dt(i));
    e(i)=max(abs(A-B));
end
% plotxx(x,e,'semilogy',Dt,e,'loglog','plot')

loglog(Dt,e,'Marker','o','Color','Red','Linestyle','none')
xlabel('Time step {\delta} t'); ylabel('error')
grid on; grid minor; ax = gca;
ax.MinorGridLineStyle = '-'; ax.GridAlpha = 0.5;

fit=polyfit(log(Dt(2:m)),log(e(2:m)),1); slope = round(fit(1));

% Dx=zeros(1,m); Dt2=0.1; Dx(1)=0.1; m=7; ex=zeros(1,m); x=m:1; 
% for i=2:m
%     Dx(i)=Dx(i-1)/2;
%     A=comp(Dt2,Dx(i-1));B=comp(Dt2,Dx(i));
%     ex(i)=max(abs(A-B));
% end
% figure
% loglog(Dx,ex)
% %semilogy(x,e,'Markersize',2)
% grid on