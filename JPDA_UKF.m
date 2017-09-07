% implementing JPDA
clc
clear
close all

%% setting up the models
dt=2*60;
Tvec=0:dt:24*60*60;
NT=length(Tvec);

No=2;
f={@(tk,xk)processmodel_2body(dt,tk,xk),@(tk,xk)processmodel_2body(dt,tk,xk)};
fn=[6,6];

Ns=1;
h=@(x)measurementmodel_2body(x) ;
senspos=[0,0,0];
sensdir=[ 0.8607    0.5089    0.0156];

SIM.use='quad';
SIM.F=@(x)processmodel_jac(dt,x);
SIM.H=@(x)measurementmodel_jac(dt,x);

%% Setting JPDA properties
JPDAprops.PD=0.8;
JPDAprops.PG=0.99;
JPDAprops.Gamma=4^2; % 4 sigma
JPDAprops.lambda=1e-5;
JPDAprops.V=1000;

%% setting up the filters
xf=cell(NT,No);
Pf=cell(NT,No);
xtruth=cell(NT,No);

Xsat0=getorbits(No);
xtruth{1,1}=Xsat0(1,:)';
xtruth{1,2}=Xsat0(2,:)';

% xf{1,1}=xtruth{1,1};
% xf{1,2}=xtruth{1,2};

% xf{1,1}(1)=xf{1,1}(1)-0.5;
% xf{1,1}(3)=xf{1,1}(3)-0.05;
%
% xf{1,2}(1)=xf{1,2}(1)+0.5;
% xf{1,2}(3)=xf{1,2}(3)+0.05;

P0=blkdiag(0.01,0.01,0.01,1e-8,1e-8,1e-8);
xf{1,1}=mvnrnd(xtruth{1,1},P0)';
xf{1,2}=mvnrnd(xtruth{1,2},P0)';
% xf{1,1}=[4.0045
%     9.6177
%     0.3788
%     -0.7301];
% xf{1,2}=[   13.3862
%     10.7268
%     -0.2199
%     -0.7618];

Pf{1,1}=P0;
Pf{1,2}=P0;

hn=3;
% R=diag([(0.1)^2]);
R=diag([(0.5*pi/180)^2,(0.5*pi/180)^2,(1)^2]);
% R=diag([(0.5*pi/180)^2,(0.5*pi/180)^2]);
% R=diag([4^2,(10*pi/180)^2]);
Q={0*diag([0.01,0.01,0.01,1e-8,1e-8,1e-8]),0*diag([0.01,0.01,0.01,1e-8,1e-8,1e-8])};

a=[5,-40];
b=[20,10];
Nc=10;
% clutter = repmat(a,Nc,1) + repmat((b-a),Nc,1).*rand(Nc,2);
% data=load('Data/sim4');
% Yhist=data.JPDAprops.Yhist;
% clutter=data.clutter;

clutter_r= [5397        3191          98];
clutter=[];
while(1)
   p=mvnrnd(clutter_r,diag([1000,1000,1000])) ;
   targdir=p/norm(p);
    if acos(sensdir*targdir(:))<=60*pi/180
       clutter=[clutter;p];
    end
    if size(clutter,1)==Nc+5
        break
    end
end


%% getting the truth
C={'ro','bo'};
% figure(5)
    for k=2:1:NT
        for i=1:No
            xtruth{k,i}=f{i}(Tvec(k-1),xtruth{k-1,i});
%             plot3(xtruth{k,i}(1),xtruth{k,i}(2),xtruth{k,i}(3),C{i},'MarkerSize',10,'linewidth',2)
%             hold on
        end
%         axis([-7100,7100,-7100,7100,-7100,7100])
%         pause(1)
%         hold off
        
    end
    

% keyboard


JPDAprops.Yhist={};

%% running the filters
for k=2:NT
    disp(strcat('iter k =',num2str(k)) );
    tic
    [xf,Pf]=propagate_JPDA(xf,Pf,SIM,Tvec,k-1,k,Q,No,f,fn,'cut8');
    toc
    
    ymset={};
    for i=1:No 
        targdir=xtruth{k,1}(1:3)/norm(xtruth{k,1}(1:3));
        if acos(sensdir*targdir(:))<=60*pi/180
           ymset{i}=h(xtruth{k,i})+sqrtm(R)*randn(hn,1);
        end
        
    end
    L=length(ymset);
    if L>0
        for pp=1:Nc+5
            if L==2
                m=randi(No,1);
                ymset{L+pp}=h(xtruth{k,m})+3*sqrtm(R)*randn(hn,1);
            else
                ymset{L+pp}=h(clutter(pp,:)')+sqrtm(R)*randn(hn,1);
            end

            if length(ymset)==12
               break 
            end
        end
        tic
        [xf,Pf,JPDAprops]=MeasurementUpdate_JPDA(xf,Pf,SIM,ymset,k,R,No,h,hn,JPDAprops,'cut8');
        toc
    else
        JPDAprops.Betas{k}=cell(1,No);
        for ii=1:No
           JPDAprops.Betas{k}{ii}=zeros(13,1); 
        end
    end
%     ymset=Yhist{k};
    
    
%     figure(1)
%     plot_JPDA(xf,Pf,clutter,No,xtruth,senspos,1,k,{'r','b'},{'ro-','bo-'})
    
    %    hold on
    %    ymset{1}
    %    for pp=2:k
    %       plot(ymset{1}(1),ymset{1}(2),'r*')
    %       plot(ymset{2}(1),ymset{2}(2),'b*')
    %    end
    %    hold off
    
    %    keyboard
    
%     pause(0.2)
%     hold off
    
end

%% plotting the simulaion probabilities
BB=cell(1,No);
XF=cell(1,No);
PF=cell(1,No);
Xtruth=cell(1,No);
for i=1:No
    B=[];
    X=[];
    P=[];
    Xt=[];
    for k=2:1:NT
        B=[B;JPDAprops.Betas{k}{i}'];
        X=[X;xf{k,i}'-xtruth{k,i}'];
        P=[P;reshape(Pf{k,i},1,fn(i)^2)];
        Xt=[Xt;xtruth{k,i}'];
    end
    XF{i}=X;
    BB{i}=B;
    PF{i}=P;
    Xtruth{i}=Xt;
end

%%

figure
plot(Tvec(2:end),BB{1},'--','linewidth',2)
legend('M1^*','M2^*','M3','M4','M5','M6','M7','M8','M9','M10','M11','M12','NoMeas')
axis([Tvec(2),Tvec(end),-0.1,1.1])
xlabel('time')
ylabel('\beta(k)')
title('Association probabilities for target 1')
plot_prop_paper

figure
plot(Tvec(2:end),BB{2},'--','linewidth',2)
legend('M1^*','M2^*','M3','M4','M5','M6','M7','M8','M9','M10','M11','M12','NoMeas')
axis([Tvec(2),Tvec(end),-0.1,1.1])
xlabel('time')
ylabel('\beta(k)')
title('Association probabilities for target 2')
plot_prop_paper



figure
plot(Tvec(2:end),sqrt(sum(XF{1}.^2,2)))
xlabel('time')
ylabel('Error in states')
title('Estimation error target 1')
plot_prop_paper

figure
plot(Tvec(2:end),sqrt(sum(XF{2}.^2,2)))
xlabel('time')
ylabel('Error in states')
title('Estimation error target 2')
plot_prop_paper

figure
plot(Tvec(2:end),sqrt(sum(PF{1}.^2,2)))
xlabel('time')
ylabel('Frobenius norm of state covariance')
title('Frobenius norm for target 1')
plot_prop_paper

figure
plot(Tvec(2:end),sqrt(sum(PF{2}.^2,2)))
xlabel('time')
ylabel('Frobenius norm of state covariance')
title('Frobenius norm for target 2')
plot_prop_paper


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