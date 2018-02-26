function [xf,Pf,JPDAprops,metrics]=MeasurementUpdate_JPDA(xf,Pf,model,Radmodel,JPDAprops,metrics,ymset,k,method)



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

for radn = 1:Radmodel.Nrad
    Nm=length(ymset{radn});
    if Nm==0
        continue
    end
    
    % doing some pre computations
    PZ=cell(1,model.No);
    MZ=cell(1,model.No);
    KK=cell(1,model.No);
    
    
    if strcmp(method,'ekf')==0
        for i=1:model.No
            muprior=xf{i};
            Pprior=Pf{i}+model.artQ{i};
            [x,w]=qd_pts(muprior,Pprior);
            z=zeros(length(w),Radmodel.hn(radn));
            flg=0;
            for q=1:length(w)
                [z(q,:),gg]=Radmodel.h{radn}(x(q,:)');
                if any(isnan(gg))
                    flg=1;
                end
            end
            
            [mz,Pz]=MeanCov(z,w);
            Pz=Pz+Radmodel.R{radn};
            MZ{i}=mz;
            Pcc=CrossCov(x,muprior,z,mz,w);
            if flg==1
                PZ{i}=1e-14*eye(size(Pz));
                KK{i}=zeros(model.fn(i),Radmodel.hn(radn));
            else
                PZ{i}=Pz;
                KK{i}=Pcc/Pz;
            end
        end
        
    else
        for i=1:model.No
            muprior=xf{i};
            Pprior=Pf{i};
            H=Radmodel.h_jac{radn}(muprior);
            [mz,gg]=Radmodel.h{radn}(muprior);
            Pz=H*Pprior*H'+Radmodel.R{radn};
            MZ{i}=mz;
            
            if any(isnan(gg))
                PZ{i}=1e-14*eye(size(Pz));
                KK{i}=zeros(model.fn(i),Radmodel.hn(radn));
            else
                PZ{i}=Pz;
                KK{i}=Pprior*H'*inv(Pz);
            end
        end
    end
    
    %% JPDA Association probabilities
    % first get all events
    % Betas (i,j): i is target, j is measurement. if there are M measurements
    % then j=M+1 is the null probability
    pp = tic;
    Betas=get_JPDA_betas_fullmarginal(MZ,PZ,model,JPDAprops,ymset{radn});
    metrics.jpda_time(k) = toc(pp);
    
    Nm=length(ymset{radn});
    JPDAprops.Betas{k}=Betas;
    
    %%   Do measurement update
    for i=1:model.No
        muprior=xf{i};
        Pprior=Pf{i};
        
        v=cell(1,Nm);
        for j=1:Nm
            v{j}=ymset{radn}{j}-MZ{i};
        end
        
        Betas(i,:)
        
        inovcov=0;
        vs=0;
        for j=1:Nm
            vs=vs+Betas(i,j)*v{j};
            inovcov=inovcov+Betas(i,j)*v{j}*v{j}';
        end
        inovcov=inovcov-vs*vs';
        Ptilde=KK{i}*(inovcov)*KK{i}';
        
        Pupdated=Pprior-KK{i}*PZ{i}*KK{i}';
        
        
        xf{i}=muprior+KK{i}*vs;
        Pf{i}=Betas(i,end)*Pprior+(1-Betas(i,end))*Pupdated+Ptilde;
    end
    
end
end