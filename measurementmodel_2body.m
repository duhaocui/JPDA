function output=measurementmodel_2body(x)

[azimuth,elevation,r] = cart2sph(x(1),x(2),x(3));
output=[azimuth,elevation,r]';
% output=[azimuth,elevation]';
end