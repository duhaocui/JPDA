function [v,R,p]=vec_radar_coordchange(v,RadPosPolar,meth)
% Nrad=size(RadPos,1);
% for i=1:1:Nrad

th=RadPosPolar(1);% perp to equator
phi=RadPosPolar(2);% along equator
r=RadPosPolar(3);
x=r*cos(th)*cos(phi);
y=r*cos(th)*sin(phi);
z=r*sin(th);
p=[x;y;z];

zv=[x,y,z]/norm([x,y,z]);
yv=-cross(zv,[0,0,1]);
yv=yv/norm(yv);
xv=-cross(zv,yv);
xv=xv/norm(xv);
% keyboard
R=[xv;yv;zv];

if strcmp(meth,'ecef2local')
    v=R*(v(:)-[x;y;z]);
else
    v=R'*v(:)+[x;y;z];
    
end






% end
