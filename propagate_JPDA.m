function [xf1,Pf1]=propagate_JPDA(xf,Pf,Tvec,k,k1,model,method)
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
        qd_pts=@get_cut8_points;
    case 'gh'
        qd_pts=@(m,P)GH_pts(m,P,para);
    case 'ekf'
        qd_pts=NaN;
    otherwise
        error('smthg is wrong: DONT ask me what')
end
%     mu_ut
%     diag(P_ut)

S=cell(1,model.No);
parfor i=1:model.No
    if strcmp(method,'ekf')==0
        try
        [x,w]=qd_pts(xf{i},Pf{i});
        catch
            keyboard
        end
        Y=zeros(size(x));
        for j=1:1:length(w)
            Y(j,:)=model.f{i}(Tvec(k),x(j,:)');
        end
        
        [N,n]=size(Y);
        W=repmat(w,1,n);
        mk=sum(W.*Y,1)';
        
        
        MU=repmat(mk',N,1);
        X=Y-MU;
        Pk=X'*(W.*X);
    elseif strcmp(method,'ekf')
        F=model.f_jac(Tvec(k),xf{i});
        mk=model.f{i}(Tvec(k),xf{i});
        Pk=F*Pf{i}*F';
        
    end
    S{i}={mk,Pk+model.Q{i}};
%     xf{k1,i}=mk;
%     Pf{k1,i}=Pk+Q{i};
end

xf1=cell(size(xf));
Pf1=cell(size(Pf));

for i=1:model.No
    xf1{i}=S{i}{1};
    Pf1{i}=S{i}{2};
end