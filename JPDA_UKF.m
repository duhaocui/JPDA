% implementing JPDA

clc
clear
close all

%% constants
Re=6378.1;

%% setting up the models
dt=2*60;
Tvec=0:dt:24*60*60;
NT=length(Tvec);



Ns=5;
h=@(x)measurementmodel_2body(x) ;
senspos=[0,0,0];
sensdir=[ 0.8607    0.5089    0.0156];

model.use_quad=true;

% SIM.use='quad';
% SIM.F=@(x)processmodel_jac(dt,x);
% SIM.H=@(x)measurementmodel_jac(dt,x);

%% set cataloged objects
model.No=20;
model.Xsat0=getorbits(model.No);
model.xtruth=cell(1,model.No);

for i=1:model.No
    model.xtruth{i}=zeros(NT,6);
    model.xtruth{i}(1,:)=model.Xsat0(i,:);
    model.f{i}=@(tk,xk)processmodel_2body(dt,tk,xk);
    model.fn(i)=6;
    model.Q{i}=1e-2*diag([1e-9,1e-9,1e-9,1e-12,1e-12,1e-12]);
end


xtruth=model.xtruth;
f=model.f;
parfor i=1:model.No
    i
    for k=2:1:NT   
        xtruth{i}(k,:)=f{i}(Tvec(k-1),xtruth{i}(k-1,:));
    end
end
model.xtruth=xtruth;


%% set Uncataloged objects
model.noncatlog.No=2;
model.noncatlog.Xsat0=getorbits(model.noncatlog.No);
model.noncatlog.xtruth=cell(1,model.noncatlog.No);

for i=1:model.noncatlog.No
    model.noncatlog.xtruth{i}=zeros(NT,6);
    model.noncatlog.xtruth{i}(1,:)=model.noncatlog.Xsat0(i,:);
    model.noncatlog.f{i}=@(tk,xk)processmodel_2body(dt,tk,xk);
    model.noncatlog.fn(i)=6;
    model.noncatlog.Q{i}=1e-2*diag([1e-9,1e-9,1e-9,1e-12,1e-12,1e-12]);
end


xtruth=model.noncatlog.xtruth;
f=model.f;
parfor i=1:model.noncatlog.No
    i
    for k=2:1:NT
        
        xtruth{i}(k,:)=f{i}(Tvec(k-1),xtruth{i}(k-1,:));
    end
end
model.noncatlog.xtruth=xtruth;


%% radars (lat,long, altitude) in geod+edic
Nrad=2;

hn=2;
R={diag([(0.5*pi/180)^2,(0.5*pi/180)^2]),diag([(0.5*pi/180)^2,(0.5*pi/180)^2])};

close all
th_rng=[-10*pi/180,10*pi/180];
phi_rng=[0,pi];
% th=RadPosPolar(1);% perp to equator
% phi=RadPosPolar(2);% along equator
% r=RadPosPolar(3);

[Ang,~]=uniform_sigma_pts([th_rng(1),phi_rng(1)],[th_rng(2),phi_rng(2)],4);


PolarPositions=cell(1,Nrad);
k=1;
for i=1:1:size(Ang,1)
    PolarPositions{i}=[Ang(i,1),Ang(i,2),Re];  % [th,phi,Re] Re means on the surface
    
    k=k+1;
    if k>Nrad
        break
    end
end

SensGeom={[pi/4,10000],[pi/4,10000]};  %misc paras such as [cone_angle,max dist of meas]

Radmodel.SensGeom=SensGeom;
Radmodel.PolarPositions=PolarPositions;
Radmodel.Nrad=Nrad;
Radmodel.hn=[hn,hn];
Radmodel.h={ @(x)measurementmodel_2body(x,Radmodel.PolarPositions{1},Radmodel.SensGeom{1}),...
    @(x)measurementmodel_2body(x,Radmodel.PolarPositions{2},Radmodel.SensGeom{2}) };

Radmodel.Nclutter=10;

X1=generate_rnd_clutter(10,Radmodel.PolarPositions{1},Radmodel.SensGeom{1})
X2=generate_rnd_clutter(10,Radmodel.PolarPositions{2},Radmodel.SensGeom{2})


figure
plot_sat_radar_system_jpda(Radmodel)
hold on
for i=1:model.No
    plot3(model.xtruth{i}(:,1),model.xtruth{i}(:,2),model.xtruth{i}(:,3),'k--')
end
plot3(X1(:,1),X1(:,2),X1(:,3),'k+','linewidth',2,'MarkerSize',6)
plot3(X2(:,1),X2(:,2),X2(:,3),'k+','linewidth',2,'MarkerSize',6)



% keyboard



figure
xlabel('x')
ylabel('y')
zlabel('z')
plot_sat_radar_system_jpda(Radmodel)
hold on
% for i=1:model.No
%     plot3(model.xtruth{i}(:,1),model.xtruth{i}(:,2),model.xtruth{i}(:,3),'k--')
% end
dothdl=cell(1,model.No);
for i=1:model.No
    dothdl{i}=plot3(model.xtruth{i}(k,1),model.xtruth{i}(k,2),model.xtruth{i}(k,3),'ko','MarkerSize',6,'linewidth',2);
end


for k=1:1:NT

    k



    for i=1:model.No
        h=0;
        for j=1:Radmodel.Nrad
            hh=measurementmodel_2body(model.xtruth{i}(k,1:3)',Radmodel.PolarPositions{j},Radmodel.SensGeom{j});
            if any(isnan(hh))
                h=0;
            else
                h=1;
                break
            end
        end
        if h==0
            set(dothdl{i},'XData', model.xtruth{i}(k,1),'YData', model.xtruth{i}(k,2),'ZData', model.xtruth{i}(k,3),'Color','k');
        else
            set(dothdl{i},'XData', model.xtruth{i}(k,1),'YData', model.xtruth{i}(k,2),'ZData', model.xtruth{i}(k,3),'Color','r');
        end
    end
    

    for i=1:Radmodel.Nrad
        
        
    end
    
    pause(0.5)
    
end




%% Now remove the sattellites that are not observable
Satobserve=zeros(1,model.No);
for k=1:1:NT
    k
    h=0;
    for i=1:model.No
        for j=1:Radmodel.Nrad
            hh=measurementmodel_2body(model.xtruth{i}(k,1:3)',Radmodel.PolarPositions{j},Radmodel.SensGeom{j});
            if any(isnan(hh))
                h=0;
            else
                h=1;
                break
            end
        end
        if h==1
            Satobserve(i)=Satobserve(i)+1;
        end
    end
end

Satobserve


%% Setting JPDA properties
JPDAprops.PD=0.8;
JPDAprops.PG=0.99;
JPDAprops.Gamma=4^2; % 4 sigma
JPDAprops.lambda=1e-5;
JPDAprops.V=1000;

JPDAprops_ukf=JPDAprops;
JPDAprops_cut4=JPDAprops;
JPDAprops_cut6=JPDAprops;

%% setting up the filters
% xf=cell(NT,No);
% Pf=cell(NT,No);
% 
% P0=blkdiag(0.01,0.01,0.01,1e-8,1e-8,1e-8);
% xtruth=cell(NT,No);
% xtruth{1,1}=Xsat0(1,:)';
% xtruth{1,2}=Xsat0(2,:)';
% 
% 
% a=[5,-40];
% b=[20,10];
% Nc=10;
% % clutter = repmat(a,Nc,1) + repmat((b-a),Nc,1).*rand(Nc,2);
% % data=load('Data/sim4');
% % Yhist=data.JPDAprops.Yhist;
% % clutter=data.clutter;
% 
% clutter_r= [5397        3191          98];
% clutter=[];
% while(1)
%    p=mvnrnd(clutter_r,diag([1000,1000,1000])) ;
%    targdir=p/norm(p);
%     if acos(sensdir*targdir(:))<=60*pi/180
%        clutter=[clutter;p];
%     end
%     if size(clutter,1)==Nc+5
%         break
%     end
% end


%% getting the truth
P0=blkdiag(0.01,0.01,0.01,1e-8,1e-8,1e-8);
% xf{k,i} i=1,2,...,No
xfukf=cell(NT,model.No);
Pfukf=cell(NT,model.No);

xfcut4=cell(NT,model.No);
Pfcut4=cell(NT,model.No);

xfcut6=cell(NT,model.No);
Pfcut6=cell(NT,model.No);

for i=1:model.No
xfukf{1,i}=mvnrnd(model.xtruth{1,1},P0)';
Pfukf{1,i}=P0;

xfcut4{1,i}=mvnrnd(model.xtruth{1,1},P0)';
Pfcut4{1,i}=P0;

xfcut6{1,i}=mvnrnd(model.xtruth{1,1},P0)';
Pfcut6{1,i}=P0;
end

metrics_ukf.propagate_time=zeros(1,NT);
metrics_cut4.propagate_time=zeros(1,NT);
metrics_cut6.propagate_time=zeros(1,NT);

metrics_ukf.measurement_time=zeros(1,NT);
metrics_cut4.measurement_time=zeros(1,NT);
metrics_cut6.measurement_time=zeros(1,NT);


%% running the filters
for k=2:NT
    disp(strcat('iter k =',num2str(k)) );
    metrics_ukf.propagate_time(k)=tic;
    [xfukf,Pfukf]=propagate_JPDA(xfukf,Pfukf,Tvec,k-1,k,model,'ukf');
    metrics_ukf.propagate_time(k)=toc(metrics_ukf.propagate_time(k));
    
    metrics_cut4.propagate_time(k)=tic;
    [xfcut4,Pfcut4]=propagate_JPDA(xfcut4,Pfcut4,Tvec,k-1,k,model,'cut4');
    metrics_cut4.propagate_time(k)=toc(metrics_cut4.propagate_time(k));
    
    metrics_cut6.propagate_time(k)=tic;
    [xfcut6,Pfcut6]=propagate_JPDA(xfcut6,Pfcut6,Tvec,k-1,k,model,'cut6');
    metrics_cut6.propagate_time(k)=toc(metrics_cut6.propagate_time(k));
    
    [metrics_ukf.propagate_time(k),metrics_cut4.propagate_time(k),metrics_cut6.propagate_time(k)]
    
    ymset=cell(1,Radmodel.Nrad);
 
    for i=1:model.No
        x=model.xtruth{i}(k,1:3);
        for j=1:Radmodel.Nrad
            hh=Radmodel.h{j}(x);
            if any(isnan(hh))
               
            else
                ymset{j}{end+1}=hh;
            end
        end
    end
    % now add the clutter 
    for j=1:Radmodel.Nrad
        Xclutter=generate_rnd_clutter(Radmodel.Nclutter,Radmodel.PolarPositions{j},Radmodel.SensGeom{j});
        for s=1:size(Xclutter,1)
            if any(isnan(hh))==0
                hh=Radmodel.h{j}(Xclutter(s,:));
                ymset{j}{end+1}=hh;
            end
        end
    end
    
    L=length(ymset);
    if L>0
        metrics_ukf.measurement_time(k)=tic;
        [xfukf,Pfukf,JPDAprops_k_ukf]=MeasurementUpdate_JPDA(xfukf,Pfukf,model,Radmodel,ymset,k,'ukf');
        metrics_ukf.measurement_time(k)=toc(metrics_ukf.measurement_time(k));
        
        
    else
        JPDAprops_ukf.Betas{k}=cell(1,No);
        JPDAprops_cut4.Betas{k}=cell(1,No);
        JPDAprops_cut6.Betas{k}=cell(1,No);
        for ii=1:model.No
           JPDAprops_ukf.Betas{k}{ii}=zeros(length(ymset),1); 
           JPDAprops_cut4.Betas{k}{ii}=zeros(length(ymset),1); 
           JPDAprops_cut6.Betas{k}{ii}=zeros(length(ymset),1); 
           
        end
    end

    
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

%

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
plot(Tvec(2:end),sqrt(sum(XF{1}.^2,2)),'linewidth',2)
xlabel('time')
ylabel('Error in states')
title('Estimation error target 1')
plot_prop_paper

figure
plot(Tvec(2:end),sqrt(sum(XF{2}.^2,2)),'linewidth',2)
set(gca,'YTick',1:15)
xlabel('time')
ylabel('Error in states')
title('Estimation error target 2')
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