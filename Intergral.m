function [norm_P,I] = Intergral(alpha,w,rho,thickness,omega,c)
% Contribution of a mode to the pressure on the groud surface
% Intergral of energy flux \int P^2/rho dz / c
%% Input
% alpha: the acoustic-wave velocity in each layer, from the ground surface
%        to the top halfspace;
% w:     in-line wind speed in each layer;
% rho:   air density in each layer;
% thickness: thickness of each alyer;
% omega: angular frequency
% c:     the phase velocity of the model at the frequency
%% Output
% norm_P: the contribution of the mode
% I:     total intergral of the eigenfunction of the mode
altitude   = sum(thickness);
k = omega/c;
nlay = length(alpha)-1;
%% Halfspace
% eigenfunction at the depth for the halfspace
fh = zeros(2,1);
v  = sqrt( (1-w(end)^2/alpha(end)^2)*k^2-omega^2/alpha(end)^2 +2*w(end)*k*omega/alpha(end)^2 );
fh(1) = exp(-v*altitude);
fh(2) = -1i/(omega-w(end)*k)/rho(end)*v*exp(-v*altitude);
I = (altitude / (2*v) * fh(1)*conj(fh(1)));
if  ~isreal(v)
    fprintf('%.1f Imaginary part \n',omega/2/pi);
end
%% From the halfspace to first layer
f = fh;
% Equation 4 in Supplementary
for ilay = nlay:-1:1
    a = f(1);
    P = zeros(2,2);
    dz = -thickness(ilay);
    v  = sqrt( (1-w(ilay)^2/alpha(ilay)^2)*k^2-omega^2/alpha(ilay)^2+2*w(ilay)*k*omega/alpha(ilay)^2 );
    tm = 1i*v/omega/rho(ilay);
    P(1,1) = cosh(v*dz);
    P(1,2) = sinh(v*dz)/tm;
    P(2,1) = sinh(v*dz)*tm;
    P(2,2) = P(1,1);
    f = P*f;
    b = f(1);
    I = I + (a*conj(a)+b*conj(b))/2*abs(dz);
end
norm_P =  real(f(1)) ^2  / I;
end