function [h,G]=measurementmodel_2body(x,PolarPosition,SensGeom)

x=x(1:3);
x=x(:);
[xenu,~,~]=vec_radar_coordchange(x,PolarPosition,'ecef2local');
% x=x'-ecef_ref(Srad,:);
% xenu = ecef2enu(x*1e3,ecef_ref(Srad,:)*1e3)/1000;

xenu(:)';
r=norm(xenu);
th=atan2(xenu(1),xenu(2));
phi=atan2(sqrt(xenu(1)^2+xenu(2)^2),xenu(3));

G=1;
if abs(phi)>SensGeom(1) || r > SensGeom(2)
    G=NaN;
end

% if hn==3
%     h=[r;th;phi];
% %  h=[x(1);x(2);x(3)];
% elseif hn==2
% %     h=[th;phi];
% %   h=[x(1);x(2)];
% %  h=[r;phi];
%  h=[r;th];
% elseif hn==1
%     h=r;
% end

h=[th;phi];

% % [azimuth,elevation,r] = cart2sph(x(1),x(2),x(3));
% [v,Rot,p]=vec_radar_coordchange(x,PolarPosition,'ecef2local');
% 
% % output=[azimuth,elevation,r]';
% output=[azimuth,elevation]';
end