
function X=generate_rnd_clutter(N,PolarPositions,SensGeom)
[v,Rot,p]=vec_radar_coordchange([0 0 0],PolarPositions,'local2ecef');
Rot=Rot';
coneang = SensGeom(1);
R = SensGeom(2);

X=zeros(N,3);
for i=1:N
%     r=1000+(R-1000)*rand(1);
    theta=0+(2*pi-0)*rand(1);
    phi=-coneang+(2*coneang)*rand(1);
    r=0+(R*tan(phi)-0)*rand(1);
    x=r.*cos(theta);
    y=r.*sin(theta);
    z=r./tan(phi);
    X(i,:)=Rot*[x;y;z]+ p;
    
end



end





