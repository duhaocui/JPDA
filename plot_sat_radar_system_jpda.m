function plot_sat_radar_system_jpda(Radars)
Re=6378.1;
% figure
hold on

% end
%% plot radar pos and their cones
for i=1:1:Radars.Nrad
    [v,Rot,p]=vec_radar_coordchange([0 0 0],Radars.PolarPositions{i},'local2ecef');
    Rot=Rot';
%     keyboard
    plot3(v(1),v(2),v(3),'ko', 'MarkerFaceColor','k','MarkerSize',6)
    
    coneang=Radars.SensGeom{i}(1);
    R=Radars.SensGeom{i}(2);
%     R=Radars{i}.MaxRange;
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
    
%     [theta,phi]=meshgrid(linspace(0,coneang,30),linspace(0,2*pi,30));
%     x=R/cos(coneang)*sin(theta).*cos(phi);
%     y=R/cos(coneang)*sin(theta).*sin(phi);
%     z=R/cos(coneang)*cos(theta);
%     for j=1:1:size(x,1)
%         for h=1:1:size(x,2)
%             
%             XX=Rot*[x(j,h);y(j,h);z(j,h)]+p;
%             x(j,h)=XX(1);
%             y(j,h)=XX(2);
%             z(j,h)=XX(3);
%         end
%     end
%     surf(x,y,z,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
end

[x,y,z] = sphere;
surf(Re*x,Re*y,Re*z,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
colormap  winter
alpha(0.3)

% plot the earth centric coordinate system
plot3([0,1000],[0,0],[0,0],'r','linewidth',4)
plot3([0,0],[0,1000],[0,0],'b','linewidth',4)
plot3([0,0],[0,0],[0,1000],'g','linewidth',4)

grid on;
axis square
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
plot_prop_paper

hold off








