function plot_sat_radar_system2_anime(Satellites,Radars,MeasPairs,Constants,yplottruth,k)
Re=6378.1;
% figure
hold on
%% plot true sat trajs and current sat positions
% if isempty(yplottruth)==0
for i=1:1:Constants.Nsat
    
    xx1=yplottruth{i}(k,1);
    xx2=yplottruth{i}(k,2);
    xx3=yplottruth{i}(k,3);
    
    if any(MeasPairs{k}(i,:)==1)
        plot3(xx1,xx2,xx3,'ro', 'linewidth',1,'MarkerSize',10,'MarkerFaceColor','r');
    else
        plot3(xx1,xx2,xx3,'go', 'linewidth',1,'MarkerSize',10,'MarkerFaceColor','g')
    end
    % plot3(xx1(1),xx2(1),xx3(1),'bo', 'linewidth',1.5)
end
% end
%% plot radar pos and their cones
for i=1:1:Constants.Nrad
    [v,Rot,p]=vec_radar_coordchange([0 0 0],Radars{i}.PolarPositions,'local2ecef');
    Rot=Rot';
    
    plot3(v(1),v(2),v(3),'ko', 'MarkerFaceColor','k','MarkerSize',6)
    
    if any(MeasPairs{k}(:,i)==1)==0
        continue
    end
    
    coneang=Radars{i}.ConeAngle;
    R=Radars{i}.MaxRange;
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
    
end

[x,y,z] = sphere;
surf(Re*x,Re*y,Re*z,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
colormap  winter
alpha(0.3)

% plot the earth centric coordinate system
plot3([0,1000],[0,0],[0,0],'r','linewidth',4)
plot3([0,0],[0,1000],[0,0],'b','linewidth',4)
plot3([0,0],[0,0],[0,1000],'g','linewidth',4)

xlim([-10000,10000])
ylim([-10000,10000])
zlim([-10000,10000])

grid on;
axis square
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
plot_prop_paper

hold off








