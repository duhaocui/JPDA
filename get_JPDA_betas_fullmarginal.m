function Betas=get_JPDA_betas_fullmarginal(xf,Pf,model,Radmodel,JPDAprops,ymset)
% first get all events
% Betas (i,j): i is target, j is measurement. if there are M measurements
% then j=M+1 is the null probability

Nm=length(ymset);
T=[];
II=eye(Nm);
II=[II;zeros(1,Nm)];
for i=1:model.No
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
for i=1:1:model.No
    A=cat(3,A,T(:,(i-1)*Nm+1:i*Nm));
end


TargetIndList=1:1:model.No;
Beta=zeros(model.No,Nm+1);
V=JPDAprops.V;
PD=JPDAprops.PD;
Gamma=JPDAprops.Gamma;

% doing some pre computations
if strcmp(method,'ekf')==0
    PZ=cell(1,model.No);
    MZ=cell(1,model.No);
    for i=1:No
        muprior=xf{k,i};
        Pprior=Pf{k,i};
        [x,w]=qd_pts(muprior,Pprior);
        z=zeros(length(w),hn);
        for q=1:length(w)
            z(q,:)=model.h(x(q,:)');
        end
        
        [mz,Pz]=MeanCov(z,w);
        Pz=Pz+R;
        MZ{i}=mz;
        PZ{i}=Pz;
    end
end



for i=1:model.No
    
    
    
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




end