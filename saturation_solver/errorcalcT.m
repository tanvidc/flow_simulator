%solutions obtained by halving the Dx m-times
%log-log plot, slope gives the order of consistency!
%needs transerrorFxn
%Tanvi Chheda 26 May 2016

m=6; Dx=zeros(m,1);
Dx(1)=1; e=zeros(m,1); 
for i=1:m-1
    e(i)=transerrorFxn(Dx(i));
    Dx(i+1)=Dx(i)/2;
end
e(m)=transerrorFxn(Dx(m));

loglog(Dx,e,'Marker','o','Color','Red','Linestyle','none')
xlabel('Dx'); ylabel('error')
grid on; grid minor; 

fit=polyfit(log(Dx(2:m)),log(e(2:m)),1); slope = round(fit(1));

