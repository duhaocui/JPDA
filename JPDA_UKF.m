% implementing JPDA
clc
clear
close all

%% setting up the models
dt=0.5;
Tvec=0:dt:50;
NT=length(Tvec);

No=2;
f={@(x)processmodel(dt,x),@(x)processmodel(dt,x)};
fn=[4,4];

Ns=1;
h=@(x)measurementmodel(x) ;
senspos=[0,0];

SIM.use='ekf';
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

xtruth{1,1}=[ 5,10,0.3,-0.8]';
xtruth{1,2}=[ 13,10,-0.2,-0.8]';

% xf{1,1}=xtruth{1,1};
% xf{1,2}=xtruth{1,2};

% xf{1,1}(1)=xf{1,1}(1)-0.5;
% xf{1,1}(3)=xf{1,1}(3)-0.05;
%
% xf{1,2}(1)=xf{1,2}(1)+0.5;
% xf{1,2}(3)=xf{1,2}(3)+0.05;

P0=diag([0.5,0.5,0.01,0.01]);
% xf{1,1}=mvnrnd(xtruth{1,1},P0)';
% xf{1,2}=mvnrnd(xtruth{1,2},P0)';
xf{1,1}=[4.0045
    9.6177
    0.3788
    -0.7301];
xf{1,2}=[   13.3862
    10.7268
    -0.2199
    -0.7618];

Pf{1,1}=P0;
Pf{1,2}=P0;

hn=1;
% R=diag([(0.1)^2]);
R=diag([(0.5*pi/180)^2]);
% R=diag([0.5^2,(1*pi/180)^2]);
Q={1*diag([0.001,0.001,0.00001,0.00001]),diag([0.001,0.001,0.00001,0.00001])};

a=[5,-40];
b=[20,10];
Nc=10;
clutter = repmat(a,Nc,1) + repmat((b-a),Nc,1).*rand(Nc,2);
%% getting the truth
for i=1:No
    for k=2:1:NT
        xtruth{k,i}=f{i}(xtruth{k-1,i});
    end
end

JPDAprops.Yhist={};

%% running the filters
for k=2:NT
    
    
    [xf,Pf]=propagate_JPDA(xf,Pf,SIM,k-1,k,Q,No,f,fn,'ut');
    
    ymset={h(xtruth{k,2})+sqrtm(R)*randn(hn,1),h(xtruth{k,1})+sqrtm(R)*randn(hn,1) };
    for pp=No+1:No+size(clutter,1)
        ymset{pp}=h(clutter(pp-No,:)')+sqrtm(R)*randn(hn,1);
    end
    
    [xf,Pf,JPDAprops]=MeasurementUpdate_JPDA(xf,Pf,SIM,ymset,k,R,No,h,hn,JPDAprops,'ut');
    
    figure(1)
    plot_JPDA(xf,Pf,clutter,No,xtruth,senspos,1,k,{'r','b'},{'ro-','bo-'})
    
    %    hold on
    %    ymset{1}
    %    for pp=2:k
    %       plot(ymset{1}(1),ymset{1}(2),'r*')
    %       plot(ymset{2}(1),ymset{2}(2),'b*')
    %    end
    %    hold off
    
    %    keyboard
    
    pause(0.2)
    hold off
    
end

%% plotting the simulaion probabilities
BB=cell(1,No);
XF=cell(1,No);
PF=cell(1,No);

for i=1:No
    B=[];
    X=[];
    P=[];
    for k=2:1:NT
        B=[B;JPDAprops.Betas{k}{i}'];
        X=[X,xf{k,i}'-xtruth{k,i}'];
        P=[P,reshape(Pf{k,i},1,fn^2)];
    end
    XF{i}=X;
    BB{i}=B;
    PF{i}=P;
end


figure
plot(Tvec(2:end),BB{1},'--','linewidth',2)
legend('M1^*','M2^*','M3','M4','M5','M6','M7','M8','M9','M10','M11','M12','NoMeas')
axis([Tvec(2),Tvec(end),-0.1,1.1])
plot_prop_paper

figure
plot(Tvec(2:end),BB{2},'--','linewidth',2)
legend('M1^*','M2^*','M3','M4','M5','M6','M7','M8','M9','M10','M11','M12','NoMeas')
axis([Tvec(2),Tvec(end),-0.1,1.1])
plot_prop_paper



figure
plot(Tvec(2:end),sqrt(sum(XF{1}.^2,2)))
plot_prop_paper

figure
plot(Tvec(2:end),sqrt(sum(XF{2}.^2,2)))
plot_prop_paper

figure
plot(Tvec(2:end),sqrt(sum(PF{1}.^2,2)))
plot_prop_paper

figure
plot(Tvec(2:end),sqrt(sum(PF{2}.^2,2)))
plot_prop_paper



