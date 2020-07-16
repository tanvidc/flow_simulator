% Script to find, for different viscocities of oil and water:
% 1. Position of shock at times 2. Time for shock to reach other end
% 3. Saturation at the end with time

% viscratio=[0.2, 1, 5, 25]; 
% [~, n]=size(viscratio);
% 
% shloc=zeros(5000,n);
% reachtime=zeros(n,1); 
t=1:1000;
% figure(1); hold on
% for v=1:n
%     [~, shloc(:,v)]=shkfront(viscratio(v));
% %     for i=1:5000
% %         if shloc(i,v)==999
% %             reachtime(v)=i;
% %             break 
% %         end
% %     end
% end

% %Position of shock in x and t, slope gives shock velocity.
% plot(t,shloc(:,1)/10,t,shloc(:,2)/10,t,shloc(:,3)/10,...
%     t,shloc(:,4)/10,'LineWidth',1.5)
% xlabel('Time'); ylabel('Position of shock front')
% legend('0.2','1','5','25','Location','southeast')

%Time for shock to reach other end
% loglog(viscratio,reachtime)

[s,~]=shkfront(10);
satend=s(999,:);
plot(satend,'LineWidth',2)
xlabel('Time step'); ylabel('Saturation at producer')
ylim([0,1]);