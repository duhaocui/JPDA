% implementing JPDA

function JPDA_singlesense(NN)

% clc
% clear
% close all

%% constants
Re=6378.1;

%% setting up the models
dt=2*60;
Tvec=0:dt:24*60*60;
NT=length(Tvec);

method.plot_anime=false;
MU=398601.2; %km^3/sec^2

%% set cataloged objects
model.No=20;
Xsat0=getorbits(model.No);
xtruth=cell(1,model.No);

for i=1:model.No
    xtruth{i}=zeros(NT,6);
    xtruth{i}(1,:)=Xsat0(i,:);
    model.f{i}=@(tk,xk)processmodel_2body(dt,MU,tk,xk);
    model.fn(i)=6;
    model.Q{i}=1e-10*diag([1e-4,1e-4,1e-4,1e-8,1e-8,1e-8]);
    model.artQ{i}=1e0*diag([1e-6,1e-6,1e-6,1e-9,1e-9,1e-9]);
end



f=model.f;
for i=1:model.No
    for k=2:1:NT
        xtruth{i}(k,:)=f{i}(Tvec(k-1),xtruth{i}(k-1,:));
    end
end



%% set Uncataloged objects
model.noncatlog.No=2;
noncatlog.Xsat0=getorbits(model.noncatlog.No);
noncatlog.xtruth=cell(1,model.noncatlog.No);

for i=1:model.noncatlog.No
    noncatlog.xtruth{i}=zeros(NT,6);
    noncatlog.xtruth{i}(1,:)=noncatlog.Xsat0(i,:);
    model.noncatlog.f{i}=@(tk,xk)processmodel_2body(dt,tk,xk);
    model.noncatlog.fn(i)=6;
    %     model.noncatlog.Q{i}=1e-2*diag([1e-9,1e-9,1e-9,1e-12,1e-12,1e-12]);
end


f=model.f;
for i=1:model.noncatlog.No
    i
    for k=2:1:NT
        
        noncatlog.xtruth{i}(k,:)=f{i}(Tvec(k-1),noncatlog.xtruth{i}(k-1,:));
    end
end



%% radars (lat,long, altitude) in geod+edic
Nrad=1;

hn=2;

close all
th_rng=[-10*pi/180,10*pi/180];
phi_rng=[0,pi];
% th=RadPosPolar(1);% perp to equator
% phi=RadPosPolar(2);% along equator
% r=RadPosPolar(3);

[Ang,~]=uniform_sigma_pts([th_rng(1),phi_rng(1)],[th_rng(2),phi_rng(2)],4);



k=1;
PolarPositions={[0,-pi/6,Re]};

SensGeom={[pi/4,10000]};  %misc paras such as [cone_angle,max dist of meas]

Radmodel.R={diag([(0.1*pi/180)^2,(0.1*pi/180)^2])};

Radmodel.SensGeom=SensGeom;
Radmodel.PolarPositions=PolarPositions;
Radmodel.Nrad=Nrad;
Radmodel.hn=[hn];
Radmodel.h={ @(x)measurementmodel_2body(x,Radmodel.PolarPositions{1},Radmodel.SensGeom{1}) };

Radmodel.Nclutter=5;

X1=generate_rnd_clutter(10,Radmodel.PolarPositions{1},Radmodel.SensGeom{1})


if method.plot_anime
    figure
    plot_sat_radar_system_jpda(Radmodel)
    hold on
    for i=1:model.No
        plot3(xtruth{i}(:,1),xtruth{i}(:,2),xtruth{i}(:,3),'k--')
    end
    plot3(X1(:,1),X1(:,2),X1(:,3),'k+','linewidth',2,'MarkerSize',6)
end

% keyboard


if method.plot_anime
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
        dothdl{i}=plot3(xtruth{i}(k,1),xtruth{i}(k,2),xtruth{i}(k,3),'ko','MarkerSize',6,'linewidth',2);
    end
    
    
    for k=1:1:NT
        k
        for i=1:model.No
            h=0;
            for j=1:Radmodel.Nrad
                [hh,gg]=measurementmodel_2body(xtruth{i}(k,1:3)',Radmodel.PolarPositions{j},Radmodel.SensGeom{j});
                if any(isnan(gg))
                    h=0;
                else
                    h=1;
                    break
                end
            end
            if h==0
                set(dothdl{i},'XData', xtruth{i}(k,1),'YData', xtruth{i}(k,2),'ZData', xtruth{i}(k,3),'Color','k');
            else
                set(dothdl{i},'XData', xtruth{i}(k,1),'YData', xtruth{i}(k,2),'ZData', xtruth{i}(k,3),'Color','r');
            end
        end
        
        
        for i=1:Radmodel.Nrad
            
            
        end
        
        pause(0.1)
        
    end
end



%% Now remove the sattellites that are not observable
Satobserve=zeros(1,model.No);
Satobserve_all=zeros(model.No,Radmodel.Nrad);
for k=1:1:NT
    k
    h=0;
    for i=1:model.No
        for j=1:Radmodel.Nrad
            [hh,gg]=measurementmodel_2body(xtruth{i}(k,1:3)',Radmodel.PolarPositions{j},Radmodel.SensGeom{j});
            if any(isnan(gg))
                h=0;
            else
                h=1;
                Satobserve_all(i,j)= Satobserve_all(i,j)+1;
                break
            end
        end
        if h==1
            Satobserve(i)=Satobserve(i)+1;
        end
    end
end

Satobserve

remtargs=[];
goodtargs=[];
for i=1:length(Satobserve)
    if Satobserve(i)<50
        remtargs=[remtargs,i];
    else
        goodtargs=[goodtargs,i];
    end
end
remtargs=[remtargs,goodtargs(4:end)];

model.f(remtargs)=[];
model.fn(remtargs)=[];
model.Q(remtargs)=[];
model.artQ(remtargs)=[];
model.No=length(model.f);

xtruth(remtargs)=[];
Satobserve(remtargs)=[];

% keyboard

%% Setting JPDA properties
JPDAprops.PD=0.8;
JPDAprops.PG=0.99;
JPDAprops.Gamma=4^2; % 4 sigma
JPDAprops.lambda=1e-5;
JPDAprops.V=1000;

JPDAprops_ukf=JPDAprops;
JPDAprops_cut4=JPDAprops;
JPDAprops_cut6=JPDAprops;
JPDAprops_cut8=JPDAprops;

%% getting the truth
P0=blkdiag(0.01,0.01,0.01,1e-6,1e-6,1e-6);
% xf{k,i} i=1,2,...,No
xfukf=cell(NT);
Pfukf=cell(NT);

xfcut4=cell(NT);
Pfcut4=cell(NT);

xfcut6=cell(NT);
Pfcut6=cell(NT);

xfcut8=cell(NT);
Pfcut8=cell(NT);

for i=1:model.No
    xfukf{1}{i}=mvnrnd(xtruth{i}(1,:),P0)';
    Pfukf{1}{i}=P0;
    
    xfcut4{1}{i}=mvnrnd(xtruth{i}(1,:),P0)';
    Pfcut4{1}{i}=P0;
    
    xfcut6{1}{i}=mvnrnd(xtruth{i}(1,:),P0)';
    Pfcut6{1}{i}=P0;
    
    xfcut8{1}{i}=mvnrnd(xtruth{i}(1,:),P0)';
    Pfcut8{1}{i}=P0;
    
end

metrics_ukf.propagate_time=zeros(1,NT);
metrics_cut4.propagate_time=zeros(1,NT);
metrics_cut6.propagate_time=zeros(1,NT);
metrics_cut8.propagate_time=zeros(1,NT);

metrics_ukf.measurement_time=zeros(1,NT);
metrics_cut4.measurement_time=zeros(1,NT);
metrics_cut6.measurement_time=zeros(1,NT);
metrics_cut8.measurement_time=zeros(1,NT);

metrics_ukf.jpda_time=zeros(1,NT);
metrics_cut4.jpda_time=zeros(1,NT);
metrics_cut6.jpda_time=zeros(1,NT);
metrics_cut8.jpda_time=zeros(1,NT);

% keyboard
%% running the filters
for k=2:NT
    disp(strcat('iter k =',num2str(k)) );
    pp=tic;
    [xfukf{k},Pfukf{k}]=propagate_JPDA(xfukf{k-1},Pfukf{k-1},Tvec,k-1,k,model,'ut');
    metrics_ukf.propagate_time(k)=toc(pp);
    
    pp=tic;
    [xfcut4{k},Pfcut4{k}]=propagate_JPDA(xfcut4{k-1},Pfcut4{k-1},Tvec,k-1,k,model,'cut4');
    metrics_cut4.propagate_time(k)=toc(pp);
    
    pp=tic;
    [xfcut6{k},Pfcut6{k}]=propagate_JPDA(xfcut6{k-1},Pfcut6{k-1},Tvec,k-1,k,model,'cut6');
    metrics_cut6.propagate_time(k)=toc(pp);
    
    pp=tic;
    [xfcut8{k},Pfcut8{k}]=propagate_JPDA(xfcut8{k-1},Pfcut8{k-1},Tvec,k-1,k,model,'cut8');
    metrics_cut8.propagate_time(k)=toc(pp);
    
    %     [metrics_ukf.propagate_time(k),metrics_cut4.propagate_time(k),metrics_cut6.propagate_time(k)]
    
    ymset=cell(1,Radmodel.Nrad);
    
    pp=ones(1,Radmodel.Nrad);
    for i=1:model.No
        x=xtruth{i}(k,1:3);
        for j=1:Radmodel.Nrad
            [hh,gg]=Radmodel.h{j}(x);
            if any(isnan(gg))
                
            else
                ymset{j}{pp(j)}=hh+sqrtm(Radmodel.R{j})*randn(Radmodel.hn(j),1);
                pp(j)=pp(j)+1;
            end
        end
    end
    % now add the clutter
    for i=1:model.No
        for j=1:Radmodel.Nrad
            %             x=xtruth{i}(k,1:3);
            %             Xclutter=mvnrnd(x,blkdiag(1e5,1e5,1e5),50);
            Xclutter=generate_rnd_clutter(Radmodel.Nclutter,Radmodel.PolarPositions{j},Radmodel.SensGeom{j});
            for s=1:size(Xclutter,1)
                [hh,gg]=Radmodel.h{j}(Xclutter(s,:));
                if any(isnan(gg))==0
                    ymset{j}{pp(j)}=hh;
                    pp(j)=pp(j)+1;
                end
            end
        end
    end
    
    pp=tic;
    [xfukf{k},Pfukf{k},JPDAprops_ukf,metrics_ukf]=MeasurementUpdate_JPDA(xfukf{k},Pfukf{k},model,Radmodel,JPDAprops_ukf,metrics_ukf,ymset,k,'ut');
    metrics_ukf.measurement_time(k)=toc(pp);
    
    pp=tic;
    [xfcut4{k},Pfcut4{k},JPDAprops_cut4,metrics_cut4]=MeasurementUpdate_JPDA(xfcut4{k},Pfcut4{k},model,Radmodel,JPDAprops_cut4,metrics_cut4,ymset,k,'cut4');
    metrics_cut4.measurement_time(k)=toc(pp);
    
    pp=tic;
    [xfcut6{k},Pfcut6{k},JPDAprops_cut6,metrics_cut6]=MeasurementUpdate_JPDA(xfcut6{k},Pfcut6{k},model,Radmodel,JPDAprops_cut6,metrics_cut6,ymset,k,'cut6');
    metrics_cut6.measurement_time(k)=toc(pp);
    
    pp=tic;
    [xfcut8{k},Pfcut8{k},JPDAprops_cut8,metrics_cut8]=MeasurementUpdate_JPDA(xfcut8{k},Pfcut8{k},model,Radmodel,JPDAprops_cut8,metrics_cut8,ymset,k,'cut8');
    metrics_cut8.measurement_time(k)=toc(pp);
    
    %     disp('Metrics')
    %     [metrics_ukf.propagate_time(k),metrics_cut4.propagate_time(k),metrics_cut6.propagate_time(k);
    %     metrics_ukf.jpda_time(k),metrics_cut4.jpda_time(k),metrics_cut6.jpda_time(k);
    %     metrics_ukf.measurement_time(k),metrics_cut4.measurement_time(k),metrics_cut6.measurement_time(k);   ]
    %
end

%% plotting the simulaion probabilities
XFukf=cell(1,model.No);
PFukf=cell(1,model.No);
for i=1:model.No
    B=[];
    X=[];
    P=[];
    Xt=[];
    for k=1:1:NT
        X=[X;xfukf{k}{i}(:)'-xtruth{i}(k,:)];
        P=[P;reshape(Pfukf{k}{i},1,model.fn(i)^2)];
    end
    XFukf{i}=X;
    PFukf{i}=P;
end


XFcut4=cell(1,model.No);
PFcut4=cell(1,model.No);
for i=1:model.No
    X=[];
    P=[];
    for k=1:1:NT
        X=[X;xfcut4{k}{i}(:)'-xtruth{i}(k,:)];
        P=[P;reshape(Pfcut4{k}{i},1,model.fn(i)^2)];
    end
    XFcut4{i}=X;
    PFcut4{i}=P;
end


XFcut6=cell(1,model.No);
PFcut6=cell(1,model.No);
for i=1:model.No
    X=[];
    P=[];
    for k=1:1:NT
        X=[X;xfcut6{k}{i}(:)'-xtruth{i}(k,:)];
        P=[P;reshape(Pfcut6{k}{i},1,model.fn(i)^2)];
    end
    XFcut6{i}=X;
    PFcut6{i}=P;
end

XFcut8=cell(1,model.No);
PFcut8=cell(1,model.No);
for i=1:model.No
    X=[];
    P=[];
    for k=1:1:NT
        X=[X;xfcut8{k}{i}(:)'-xtruth{i}(k,:)];
        P=[P;reshape(Pfcut8{k}{i},1,model.fn(i)^2)];
    end
    XFcut8{i}=X;
    PFcut8{i}=P;
end


Final_methods={'ukf','cut4','cut6','cut8'};
Final_solutions={{XFukf,PFukf},{XFcut4,PFcut4},{XFcut6,PFcut6},{XFcut8,PFcut8}};

%%
close all


pC={'r','b','g','k'};
pS={'o','s','*','+','^','<','>'};
sat_states={'x','y','z','vx','vy','vz'};
ERRORS=zeros(length(sat_states),length(Final_methods),model.No);
Final_methods
for i=1:model.No
    for meth=1:length(Final_methods)
        XF=Final_solutions{meth}{1}{i};
        PF=Final_solutions{meth}{2}{i};
        
        sigmas=zeros(NT,model.fn(i));
        for k=1:NT
            P=sqrtm( reshape(PF(k,:),model.fn(i),model.fn(i)));
            sigmas(k,:)=diag(P)';
        end
        for s=1:model.fn(i)
            
            %             figure(s+model.fn(max(1,i-1))*(i-1) )
            %             plot(Tvec(2:end),sqrt(sum(XF(2:end,s).^2,2)),[pS{meth},'--'],'linewidth',2)
            
            ERRORS(s,meth,i)=sqrt(sum(XF(2:end,s).^2)/length(Tvec(2:end)));
            %             hold on
            %             plot(Tvec(2:end),sqrt(sum(XF(2:end,s).^2,2))+3*sigmas(2:end,s),'--','linewidth',2)
            %             plot(Tvec(2:end),sqrt(sum(XF(2:end,s).^2,2))-3*sigmas(2:end,s),'--','linewidth',2)
            
            %             xlabel('time')
            %             ylabel(['Error in states ',sat_states{s}])
            %             title(['Satellite ',num2str(i)])
            %             legend(Final_methods{1:meth})
            %             plot_prop_paper_jpda
            
        end
    end
end

sqrt( sum(ERRORS.^2,3)/size(ERRORS,3) )


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

save(['JPDA_satsinglesense_',num2str(NN),'.mat'])