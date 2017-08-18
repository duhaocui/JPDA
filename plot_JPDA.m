function plot_JPDA(xf,Pf,clutter,No,xtruth,senspos,k0,k1,C,CF)



holdset=0;
for i=1:No
    
    XT=zeros(k1-k0,length(xtruth{k0,i}) );
    for k=k0:k1
       XT(k-k0+1,:)= xtruth{k,i};
    end
    
    plot(XT(:,1),XT(:,2),C{i},'linewidth',1)
    
    if holdset==0
        hold on;
        holdset=1;
    end
    
    XF=zeros(k1-k0,length(xf{k0,i}) );
    for k=k0:k1
       XF(k-k0+1,:)= xf{k,i};
    end
    plot(XF(:,1),XF(:,2),CF{i},'linewidth',2)
    
    plot_nsigellip(xf{k1,i}(1:2),Pf{k1,i}(1:2,1:2),2,C{i},1)
    
    
    
    
end

plot(senspos(1),senspos(2),'^','MarkerSize',15,'linewidth',2)
   plot(clutter(:,1),clutter(:,2),'k+','MarkerSize',15,'linewidth',2) 


% axis([0,20,-30,30])
hold off



end