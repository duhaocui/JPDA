mu=[10,10]
plot_nsigellip([10,10],[3,1;1,3],1,'r',1)
hold on
plot(mu(1),mu(2),'ro','MarkerSize',10,'linewidth',2)
a=[0,0];
b=[20,20];
Nc=4;
clutter = repmat(a,Nc,1) + repmat((b-a),Nc,1).*rand(Nc,2);

plot(0,0,'^','MarkerSize',15,'linewidth',2)
plot(clutter(:,1),clutter(:,2),'k+','MarkerSize',15,'linewidth',2) 
axis([-2,20,-2,20])