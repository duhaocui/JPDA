function Betas=get_JPDA_betas_fullmarginal(MZ,PZ,model,JPDAprops,ymset,out_of_gated_meas_per_target)
% first get all events
% Betas (i,j): i is target, j is measurement. if there are M measurements
% then j=M+1 is the null probability

Nm=length(ymset);
A=jpda_full_event_combos(model.No,Nm,out_of_gated_meas_per_target);
% A(i,j,k) ... k is the event count or id.
% i is the target index. and j is the measurement index
% keyboard


TargetIndList=1:1:model.No;
Betas=zeros(model.No,Nm+1);

V=JPDAprops.V;
PD=JPDAprops.PD;
Gamma=JPDAprops.Gamma;

No=model.No;
Prob_events = zeros(size(A,3),1);

for e=1:1:size(A,3) % go throught all thee events
%     if A(i,j,e)==1 % if in event e, jth meas is paired with target i
        
        phi=Nm-sum(sum(A(:,:,e),1)); % number of un associated measurements
        
        Tprod=zeros(1,No); % over targets
        Mprod=zeros(1,Nm); % over measurements
        for ei=1:1:No
            delt=any(A(ei,:,e)==1);
            Tprod(ei)=PD^delt * (1-PD)^(1-delt);
        end
        for ej=1:1:Nm  % check validation here
            tauj=any(A(:,ej,e)==1);
            if tauj==1 % if meas ej is actually used
                target_ass=TargetIndList( A(:,ej,e)==1 );
                mz=MZ{target_ass};
                Pz=PZ{target_ass};
                
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
        
        Prob_events(e)=(factorial(phi)/V^phi)*prod(Mprod)*prod(Tprod);
%         Betas(i,j)=Betas(i,j)+Pevent;

    
    
end

% keyboard

for i=1:No
    for j=1:Nm
        Betas(i,j)=sum(Prob_events(A(i,j,:)==1) ) ;
    end    
    
    % Now get the null probability
    for e=1:size(A,3)
       if sum(A(i,:,e))==0
          Betas(i,end) = Betas(i,end) + Prob_events(e);
       end
    end
    Betas(i,:) = Betas(i,:)/sum(Betas(i,:));
end




end