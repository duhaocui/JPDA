function A=jpda_full_event_combos(Nt,Nm,out_of_gated_meas_per_target)
II=eye(Nm);
II=[II;zeros(1,Nm)];
measindex=1:Nm;

%% Roughly count the number of events that have to be removed when out_of_gated_meas_per_target is given
C=0;
for i=1:Nt
    S=0;
    if length(out_of_gated_meas_per_target{i})>0
       S=count_events_fix_1(Nt-1,Nm-1);
    end
    C=C+S*length(out_of_gated_meas_per_target{i});
end

%% now find all the events
a=ones(1,Nt);
% A=zeros(Nt,Nm,min(size(II,1)-C,size(II,1))^Nt);min(size(II,1)-C,size(II,1))
A=zeros(Nt,Nm,max(min(size(II,1)-C,size(II,1)),0)^Nt);
k=1;
flg=0;
while 1
    C=zeros(Nt,Nm);
    for i=1:Nt
        C(i,:)=II(a(i),:)';
    end
    skipevent=false;
    for i=1:Nt
        if a(i)<=Nm
            if any(measindex(a(i))==out_of_gated_meas_per_target{i})
                skipevent=true;
                break;
            end 
        end
    end
    if skipevent==true
       
    else
        A(:,:,k)=C;
        k=k+1;
    end
    
    
    
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


