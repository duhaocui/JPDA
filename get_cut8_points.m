function [X,w]=get_cut8_points(mu,P)
n=length(mu);
filename=['CUT8_',num2str(n),'.mat'];
if exist(filename,'file')==2
    load(filename)
else
    [X,w]=conjugate_dir_gausspts_till_8moment(zeros(n,1),eye(n));
    save(filename,'X','w')
end

A=sqrtm(P);
for i=1:1:length(w)
X(i,:)=A*X(i,:)'+mu;
end