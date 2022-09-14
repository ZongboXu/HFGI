function w0 = PropMatW(alpha,w,rho,thickness,omega,c)
% Propagation matrix from the halfspace to the ground surface
% the vertical velocity at the surface should be 0
%% Input
% alpha: the acoustic-wave velocity in each layer, from the ground surface
%        to the top halfspace;
% w:     in-line wind speed in each layer;
% rho:   air density in each layer;
% thickness: thickness of each alyer;
% omega: angular frequency
% c:     the phase velocity of the model at the frequency
%% Output
% w0:    the vertical particle velocity on the ground surface. 
% If w0 = 0, a root is found for the frequency
altitude   = sum(thickness);
k = omega/c;
nlay = length(alpha)-1;
%% Halfspace
% eigenfunction at the depth for the halfspace
fh = zeros(2,1);
v  = sqrt( (1-w(end)^2/alpha(end)^2)*k^2-omega^2/alpha(end)^2 +2*w(end)*k*omega/alpha(end)^2 );
fh(1) = exp(-v*altitude);
fh(2) = -1i/(omega-w(end)*k)/rho(end)*v*exp(-v*altitude);
%% From the halfspace to first layer
f = fh;
for ilay = nlay:-1:1
    P = zeros(2,2);
    dz = -thickness(ilay);
    v  = sqrt( (1-w(ilay)^2/alpha(ilay)^2)*k^2-omega^2/alpha(ilay)^2+2*w(ilay)*k*omega/alpha(ilay)^2 );
    tm = 1i*v/omega/rho(ilay);
    P(1,1) = cosh(v*dz);
    P(1,2) = sinh(v*dz)/tm;
    P(2,1) = sinh(v*dz)*tm;
    P(2,2) = P(1,1);
    f = P*f;
end
% vertical velocity at the ground surface
w0 = imag(f(2));
end