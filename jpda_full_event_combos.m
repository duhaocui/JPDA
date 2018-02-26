function A=jpda_full_event_combos(Nt,Nm)
II=eye(Nm);
II=[II;zeros(1,Nm)];


a=ones(1,Nt);
A=zeros(Nt,Nm,size(II,1)^Nt);
k=1;
flg=0;
while 1
    C=zeros(Nt,Nm);
    for i=1:Nt
        C(i,:)=II(a(i),:);
    end
    A(:,:,k)=C;
    
    a(1)=a(1)+1;
    for i=1:Nt
        if a(i)>Nm+1
            if i==Nt
               flg=1; 
               break
            end
            a(i)=1;
            a(i+1)=a(i+1)+1;
        end
    end
    
    if flg==1
        break
    end
    
    k=k+1;
    
    
end
A(:,:,k+1:end)=[];

%% constraint checking
remove_index=[];
for i=1:size(A,3)
    %sum along rows
    if any(sum(A(:,:,i),2)>1) || any(sum(A(:,:,i),1)>1)
        remove_index=horzcat(remove_index,i);
    end
    
    
end
if length(remove_index)>0
    A(:,:,remove_index)=[];
end


