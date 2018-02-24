function plot_sat_radar_system_dynamic(yplottruth,Radmodel,further)
Re=6378.1;
% figure
hold on
%% plot true sat trajs and current sat positions 
if isempty(yplottruth)==0
for i=1:1:Radmodel.Nsat
    
xx1=yplottruth{i,1}(:,1);
xx2=yplottruth{i,1}(:,2);
xx3=yplottruth{i,1}(:,3);

if max(strcmpi(further,'Sathighlight'))==1
    if i==GetCellOptions({'Sathighlight',11},'Sathighlight')
        plot3(xx1,xx2,xx3,'k:', 'linewidth',1.5)
    else
        plot3(xx1,xx2,xx3,'r:', 'linewidth',1.5)
    end
end
% plot3(xx1(1),xx2(1),xx3(1),'bo', 'linewidth',1.5)
end
end
%% plot radar pos and their cones 
for i=1:1:Radmodel.Nrad
   [v,Rot,p]=vec_radar_coordchange([0 0 0],i,Radmodel.RadPos,'local2ecef');
   Rot=Rot';
   p';
plot3(v(1),v(2),v(3),'ko', 'MarkerFaceColor','k','MarkerSize',6)
% if isempty(find(MeasPairs{k}(:,2)==i))
%     continue
% end

% x1 = ecef2enu([1,0,0],Radmodel.pos(i,:)*1e3);
% x2 = ecef2enu([0,1,0],Radmodel.pos(i,:)*1e3);
% x3 = ecef2enu([0,0,1],Radmodel.pos(i,:)*1e3);
% Rot=[x1(:),x2(:),x3(:)];
% Rot=Rot';
% 
% %% have to shift and rotate the cones
% keyboard
coneang=Radmodel.SensParas(i,1);
R=Radmodel.SensParas(i,2);
r=linspace(0,R*tan(coneang),30);
theta=linspace(0,2*pi,30);
[r,theta]=meshgrid(r,theta);
x=r.*cos(theta);
y=r.*sin(theta);
z=r./tan(coneang);
for j=1:1:size(x,1)
    for h=1:1:size(x,2)
   XX=Rot*[x(j,h);y(j,h);z(j,h)]+ p;    
   x(j,h)=XX(1);
   y(j,h)=XX(2);
   z(j,h)=XX(3);
    end
end
surf(x,y,z,'EdgeColor','none','LineStyle','none','FaceLighting','phong')

[theta,phi]=meshgrid(linspace(0,coneang,30),linspace(0,2*pi,30));
x=R/cos(coneang)*sin(theta).*cos(phi);
y=R/cos(coneang)*sin(theta).*sin(phi);
z=R/cos(coneang)*cos(theta);
for j=1:1:size(x,1)
    for h=1:1:size(x,2)

   XX=Rot*[x(j,h);y(j,h);z(j,h)]+p;    
   x(j,h)=XX(1);
   y(j,h)=XX(2);
   z(j,h)=XX(3);
    end
end
% surf(x,y,z,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
% alpha(0.3)




end

[x,y,z] = sphere;
surf(Re*x,Re*y,Re*z,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
colormap  winter
alpha(0.3)




% axis(5*Re*[-1,1,-1,1,-1,1])
axis square

hold off








