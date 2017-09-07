function [xf,Pf]=propagate_JPDA(xf,Pf,SIM,Tvec,k,k1,Q,No,f,fn,method)
if k1-k~=1
    error('wrong k1 and k')
end


% global kappa
% kappa=para;
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
%     mu_ut
%     diag(P_ut)
USE=SIM.use;
S=cell(1,No);
parfor i=1:No
    if strcmp(USE,'quad')
        [x,w]=qd_pts(xf{k,i},Pf{k,i});
        
        Y=zeros(size(x));
        for j=1:1:length(w)
            Y(j,:)=f{i}(Tvec(k),x(j,:)');
        end
        
        [N,n]=size(Y);
        W=repmat(w,1,n);
        mk=sum(W.*Y,1)';
        
        
        MU=repmat(mk',N,1);
        X=Y-MU;
        Pk=X'*(W.*X);
    elseif strcmp(USE,'ekf')
        F=SIM.F(xf{k,i});
        mk=f{i}(xf{k,i});
        Pk=F*Pf{k,i}*F';
        
    end
    S{i}={mk,Pk+Q{i}};
%     xf{k1,i}=mk;
%     Pf{k1,i}=Pk+Q{i};
end

for i=1:No
xf{k1,i}=S{i}{1};
Pf{k1,i}=S{i}{2};
end