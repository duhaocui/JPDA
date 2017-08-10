function output=measurementmodel(x)

if x(1)==0
    th=0;
else
    th=atan2(x(2),x(1));
end

r=sqrt(x(1)^2+x(2)^2);

% r=x(1);
% th=x(2);

% output=[r,th]';
% output=[r]';
output=[th]';
end