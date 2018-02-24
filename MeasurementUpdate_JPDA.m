function [xf,Pf,JPDAprops]=MeasurementUpdate_JPDA(xf,Pf,model,Radmodel,ymset,k,method)



switch lower(method)
    case 'ckf'
        qd_pts=@cubature_KF_points;
    case 'ut'
        qd_pts=@(m,P)UT_sigmapoints(m,P,2);
    case 'cut4'
        qd_pts=@conjugate_dir_gausspts;
    case 'cut6'
        qd_pts=@conjugate_dir_gausspts_till_6moment_scheme2;
    case 'cut8'
        qd_pts=@conjugate_dir_gausspts_till_8moment;
    case 'gh'
        qd_pts=@(m,P)GH_pts(m,P,para);
    case 'ekf'
        qd_pts=NaN;
    otherwise
        error('smthg is wrong: DONT ask me what')
        %         qd_pts
end

%% JPDA Association probabilities
% first get all events
% Betas (i,j): i is target, j is measurement. if there are M measurements
% then j=M+1 is the null probability
Betas=get_JPDA_betas_fullmarginal(model,Radmodel);


Nm=length(ymset);

JPDAprops.Yhist{k}=ymset;
JPDAprops.pdfZ{k}=cell(1,No);
JPDAprops.Betas{k}=Betas;

% keyboard

%%   Do measurement update
for i=1:model.No
    
    if strcmp(method,'quad')==0
        muprior=xf{k,i};
        Pprior=Pf{k,i};
        [x,w]=qd_pts(muprior,Pprior);
        z=zeros(length(w),hn);
        for j=1:length(w)
            z(j,:)=model.h(x(j,:)');
        end
        
        [mz,Pz]=MeanCov(z,w);
        Pz=Pz+R;
        
        Pcc=CrossCov(x,muprior,z,mz,w);
        K=Pcc/Pz;
        
    elseif strcmp(method,'ekf')
        muprior=xf{k,i};
        Pprior=Pf{k,i};
        H=model.h_jac(muprior);
        mz=model.h(muprior);
        Pz=H*Pprior*H'+R;
        K=Pprior*H'*inv(Pz);
    end
    
    
    JPDAprops.pdfZ{k}{i}={mz,Pz};
    
    
    
    
    v=cell(1,Nm);
    for j=1:Nm
        v{j}=ymset{j}-mz;
    end
%     ss=sum(Beta(:,i))+Beta_null(i);
%     Beta(:,i)=Beta(:,i)/ss;
%     Beta_null(i)=Beta_null(i)/ss;
    
%     JPDAprops.Betas{k}=Betas;
    
    Betas(i,:)
    
    inovcov=0;
    vs=0;
    for j=1:Nm
        vs=vs+Betas(i,j)*v{j};
        inovcov=inovcov+Betas(i,j)*v{j}*v{j}';
    end
    inovcov=inovcov-vs*vs';
    Ptilde=K*(inovcov)*K';
    
    Pupdated=Pprior-K*Pz*K';
    %     keyboard
    
    
    xf{k,i}=muprior+K*vs;
    Pf{k,i}=Betas(i,end)*Pprior+(1-Betas(i,end))*Pupdated+Ptilde;
    %     i
end
end