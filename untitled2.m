% for i=50:-1:40
%    RunSIMJPDA('cut8',i)
%
% end



%%
Xcut6=cell(1,2);
Xut=cell(1,2);
Xcut8=cell(1,2);

for i=1:50
    UT=load(strcat('ut_',num2str(i)));
    p1=sum(UT.XF{1}.^2,2);
    p2=sum(UT.XF{2}.^2,2);
    if i==1
        Xut{1}=p1;
        Xut{2}=p2;
    else
        Xut{1}=Xut{1}+p1;
        Xut{2}=Xut{2}+p2;
    end
    
    CUT6=load(strcat('cut6_',num2str(i)));
    p1=sum(CUT6.XF{1}.^2,2);
    p2=sum(CUT6.XF{2}.^2,2);
    
    if i==1
        Xcut6{1}=p1;
        Xcut6{2}=p2;
    else
        Xcut6{1}=Xcut6{1}+p1;
        Xcut6{2}=Xcut6{2}+p2;
    end
    
    
    if i<=9
        CUT8=load(strcat('cut8_',num2str(i)));
        p1=sum(CUT8.XF{1}.^2,2);
        p2=sum(CUT8.XF{2}.^2,2);
        
        if i==1
            Xcut8{1}=p1;
            Xcut8{2}=p2;
        else
            Xcut8{1}=Xcut8{1}+p1;
            Xcut8{2}=Xcut8{2}+p2;
        end
    end
end


Xcut6{1}=sqrt(Xcut6{1}/50);
Xcut6{2}=sqrt(Xcut6{2}/50);

Xut{1}=sqrt(Xut{1}/50);
Xut{2}=sqrt(Xut{2}/50);

Xcut8{1}=sqrt(Xcut8{1}/9);
Xcut8{2}=sqrt(Xcut8{2}/9);

figure
plot(Tvec(2:end),Xut{1},'b--',Tvec(2:end),Xcut8{1},'k','linewidth',2)
xlabel('time')
ylabel('RMSE')
legend('UT','CUT8')
plot_prop_paper

figure
plot(Tvec(2:end),Xut{2},'b--',Tvec(2:end),Xcut8{2},'k','linewidth',2)
xlabel('time')
ylabel('RMSE')
xlabel('time')
ylabel('RMSE')
legend('UT','CUT8')
plot_prop_paper


% figure
% plot(Tvec(2:end),BB{1},'--','linewidth',2)
% legend('M1^*','M2^*','M3','M4','M5','M6','M7','M8','M9','M10','M11','M12','NoMeas')
% axis([Tvec(2),Tvec(end),-0.1,1.1])
% xlabel('time')
% ylabel('\beta(k)')
% title('Association probabilities for target 1')
% plot_prop_paper
%
% figure
% plot(Tvec(2:end),BB{2},'--','linewidth',2)
% legend('M1^*','M2^*','M3','M4','M5','M6','M7','M8','M9','M10','M11','M12','NoMeas')
% axis([Tvec(2),Tvec(end),-0.1,1.1])
% xlabel('time')
% ylabel('\beta(k)')
% title('Association probabilities for target 2')
% plot_prop_paper
%
%
%
% figure
% plot(Tvec(2:end),sqrt(sum(XF{1}.^2,2)),'linewidth',2)
% xlabel('time')
% ylabel('Error in states')
% title('Estimation error target 1')
% plot_prop_paper
%
% figure
% plot(Tvec(2:end),sqrt(sum(XF{2}.^2,2)),'linewidth',2)
% set(gca,'YTick',1:15)
% xlabel('time')
% ylabel('Error in states')
% title('Estimation error target 2')
% plot_prop_paper




%% Animation of tracking and pdfs in angle space
% data=load('Data/sim2')
% for k=2:data.NT
%     figure(1)
%     plot_JPDA(data.xf,data.Pf,data.clutter,data.No,data.xtruth,data.senspos,1,k,{'r','b'},{'ro-','bo-'})
%     pause(0.1)
%     saveas(gcf,strcat('Anime/',sprintf('sim%6.6d', k)),'png')
% end

%%
% close all
% C={'r','g'};
% yrng=linspace(-pi/2,pi/2,100);
% for k=2:1:data.NT
%     figure(2)
%     hold on
%     for i=1:data.No
%         mz=data.JPDAprops.pdfZ{k}{i}{1};
%         Pz=data.JPDAprops.pdfZ{k}{i}{2};
%         yrng=linspace(mz-3*sqrtm(Pz),mz+3*sqrtm(Pz),100);
%         pth=normpdf(yrng,mz,Pz);
%         plot(yrng,pth,C{i})
%     end
%     NM=length(data.JPDAprops.Yhist{k});
%     for j=1:NM
%         plot(data.JPDAprops.Yhist{k}{j},0,'k*')
%     end
%     axis([-pi/2,pi/2,-0.2,3])
% %     keyboard
%     pause(1)
%     clf
% end
