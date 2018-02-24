function XX = OE2XYZ(OE)
 a0=OE(1);
 e0=OE(2);
 E0=OE(3);
 w0=OE(4);
 i0=OE(5);
 Om0=OE(6);
 
 mue=398601.2;
%% function to determine position and velocity from orbital elements. 

r0 = a0*(1-e0*cos(E0));
xb = [a0*(cos(E0)-e0); a0*sqrt(1-e0^2)*sin(E0); 0]; 
xdotb = [-sqrt(mue*a0)*sin(E0)/r0; sqrt(mue*a0*(1-e0^2))*cos(E0)/r0; 0 ];

% Orbit initial orientation parameters
R = [cos(Om0), -sin(Om0), 0; sin(Om0), cos(Om0), 0; 0 0 1];
R = R*[1, 0, 0; 0, cos(i0), -sin(i0); 0 sin(i0) cos(i0)];
R = R*[cos(w0), -sin(w0), 0; sin(w0), cos(w0), 0; 0 0 1];


% initial position and velocity (inertial frame ECEF)
X = R*xb; Xdot = R*xdotb;   xb0 = [X;Xdot];

XX=zeros(6,1);
XX(1:3)=X;
XX(4:6)=Xdot;
