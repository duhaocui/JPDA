function [xf,Pf,JPDAprops]=MeasurementUpdate_JPDA(xf,Pf,SIM,ymset,k,R,No,h,hn,JPDAprops,method)



switch lower(method)
    case {'ckf'}
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
end

%% JPDA Association probabilities
% first get all events
Nm=length(ymset);
T=[];
II=eye(Nm);
II=[II;zeros(1,Nm)];
for i=1:No
    T=tense_prod_index(T,II);
end
remind=[];
for r=1:size(T,1)
    a=zeros(1,Nm);
    for i=1:No
        a=a+T(r,(i-1)*Nm+1:i*Nm);
    end
    
    if any(a>=2)
        remind=horzcat(remind,r);
    end
end

T(remind,:)=[];

A=[];
for i=1:1:No
    A=cat(3,A,T(:,(i-1)*Nm+1:i*Nm));
end


TargetIndList=1:1:No;
Beta=zeros(Nm,No);
Beta_null=zeros(1,No); % probability in which no measurements are associated
V=JPDAprops.V;
PD=JPDAprops.PD;
Gamma=JPDAprops.Gamma;

% doing some pre computations
if strcmp(SIM.use,'quad')
    PZ=cell(1,No);
    MZ=cell(1,No);
    for i=1:No
        muprior=xf{k,i};
        Pprior=Pf{k,i};
        [x,w]=qd_pts(muprior,Pprior);
        z=zeros(length(w),hn);
        for q=1:length(w)
            z(q,:)=h(x(q,:)');
        end
        
        [mz,Pz]=MeanCov(z,w);
        Pz=Pz+R;
        MZ{i}=mz;
        PZ{i}=Pz;
    end
end

USE=SIM.use;
% getting validation region
for i=1:No
    
    
    
    for j=1:Nm
        
        for e=1:1:size(A,1) % go throught all thee events
            if A(e,j,i)==1 % if in event e, jth meas is paired with target i
                
                phi=Nm-sum(sum(A(e,:,:),3)); % number of un associated measurements
                
                Tprod=zeros(1,No); % over targets
                Mprod=zeros(1,Nm); % over measurements
                for ei=1:1:No
                    delt=any(A(e,:,ei)==1);
                    Tprod(ei)=PD^delt * (1-PD)^(1-delt);
                end
                for ej=1:1:Nm  % check validation here
                    tauj=any(A(e,ej,:)==1);
                    if tauj==1 % if meas ej is actually used
                        target_ass=TargetIndList( A(e,ej,:)==1 );
                        
                        if strcmp(USE,'quad')
                            %                             muprior=xf{k,target_ass};
                            %                             Pprior=Pf{k,target_ass};
                            %                             [x,w]=qd_pts(muprior,Pprior);
                            %                             z=zeros(length(w),hn);
                            %                             for q=1:length(w)
                            %                                 z(q,:)=h(x(q,:)');
                            %                             end
                            %
                            %                             [mz,Pz]=MeanCov(z,w);
                            %                             Pz=Pz+R;
                            mz=MZ{target_ass};
                            Pz=PZ{target_ass};
                        elseif strcmp(USE,'ekf')
                            muprior=xf{k,target_ass};
                            Pprior=Pf{k,target_ass};
                            H=SIM.H(muprior);
                            mz=h(muprior);
                            Pz=H*Pprior*H'+R;
                        end
                        
                        if (ymset{ej}-mz)'*inv(Pz)*(ymset{ej}-mz)>=Gamma
                            Mprod(ej)=0; % outside the sigma ellipsoid, so its probability is 0, hence tis event has 0 prob
                        else
                            %                             keyboard
                            likelihood=mvnpdf(ymset{ej},mz,Pz);
                            Mprod(ej)=likelihood^tauj;
                        end
                    else
                        Mprod(ej)=1;
                    end
                    
                end
                
                Pevent=(factorial(phi)/V^phi)*prod(Mprod)*prod(Tprod);
                Beta(j,i)=Beta(j,i)+Pevent;
            end
            
            
        end
    end
    
    % getting the null marginal i.e. marginal events where no measurements are associated with target i
    for e=1:1:size(A,1) % go throught all thee events
        if sum(A(e,:,i))==0 % events where no meas assigned to target i
            phi=Nm-sum(sum(A(e,:,:),3)); % number of un associated measurements
            
            Tprod=zeros(1,No); % over targets
            Mprod=zeros(1,Nm); % over measurements
            for ei=1:1:No
                delt=any(A(e,:,ei)==1);
                Tprod(ei)=PD^delt * (1-PD)^(1-delt);
            end
            for ej=1:1:Nm  % check validation here
                tauj=any(A(e,ej,:)==1);
                if tauj==1 % if meas ej is actually used
                    target_ass=TargetIndList( A(e,ej,:)==1 );
                    
                    if strcmp(SIM.use,'quad')
                        %                         muprior=xf{k,target_ass};
                        %                         Pprior=Pf{k,target_ass};
                        %                         [x,w]=qd_pts(muprior,Pprior);
                        %                         z=zeros(length(w),hn);
                        %                         for q=1:length(w)
                        %                             z(q,:)=h(x(q,:)');
                        %                         end
                        %
                        %                         [mz,Pz]=MeanCov(z,w);
                        %                         Pz=Pz+R;
                        mz=MZ{target_ass};
                        Pz=PZ{target_ass};
                    elseif strcmp(SIM.use,'ekf')
                        muprior=xf{k,target_ass};
                        Pprior=Pf{k,target_ass};
                        H=SIM.H(muprior);
                        mz=h(muprior);
                        Pz=H*Pprior*H'+R;
                    end
                    
                    if (ymset{ej}-mz)'*inv(Pz)*(ymset{ej}-mz)>=Gamma
                        Mprod(ej)=0; % outside the sigma ellipsoid, so its probability is 0, hence tis event has 0 prob
                    else
                        likelihood=mvnpdf(ymset{ej},mz,Pz);
                        Mprod(ej)=likelihood^tauj;
                    end
                else
                    Mprod(ej)=1;
                end
                
            end
            
            Pevent=(factorial(phi)/V^phi)*prod(Mprod)*prod(Tprod);
            Beta_null(i)=Beta_null(i)+Pevent;
            
            
        end
    end
    
    
    
    
end

JPDAprops.Yhist{k}=ymset;
JPDAprops.pdfZ{k}=cell(1,No);
JPDAprops.Betas{k}=cell(1,No);

% keyboard

%%   Do measurement update
for i=1:No
    
    if strcmp(SIM.use,'quad')
        muprior=xf{k,i};
        Pprior=Pf{k,i};
        [x,w]=qd_pts(muprior,Pprior);
        z=zeros(length(w),hn);
        for j=1:length(w)
            z(j,:)=h(x(j,:)');
        end
        
        [mz,Pz]=MeanCov(z,w);
        Pz=Pz+R;
        
        Pcc=CrossCov(x,muprior,z,mz,w);
        K=Pcc/Pz;
        
    elseif strcmp(SIM.use,'ekf')
        muprior=xf{k,i};
        Pprior=Pf{k,i};
        H=SIM.H(muprior);
        mz=h(muprior);
        Pz=H*Pprior*H'+R;
        K=Pprior*H'*inv(Pz);
    end
    
    
    JPDAprops.pdfZ{k}{i}={mz,Pz};
    
    
    
    
    v=cell(1,Nm);
    for j=1:Nm
        v{j}=ymset{j}-mz;
    end
    ss=sum(Beta(:,i))+Beta_null(i);
    Beta(:,i)=Beta(:,i)/ss;
    Beta_null(i)=Beta_null(i)/ss;
    
    JPDAprops.Betas{k}{i}=[Beta(:,i);Beta_null(i)];
    
    Beta
    Beta_null
    
    inovcov=0;
    vs=0;
    for j=1:Nm
        vs=vs+Beta(j,i)*v{j};
        inovcov=inovcov+Beta(j,i)*v{j}*v{j}';
    end
    inovcov=inovcov-vs*vs';
    Ptilde=K*(inovcov)*K';
    
    Pupdated=Pprior-K*Pz*K';
    %     keyboard
    
    
    xf{k,i}=muprior+K*vs;
    Pf{k,i}=Beta_null(i)*Pprior+(1-Beta_null(i))*Pupdated+Ptilde;
    %     i
end
end