function [pvcomp,pvcomph] = Compliance_Mlayer( frelist, clist,Vp,Vs,den,z)
% Calculating the compliance for a 1D M-layered subsurface velocity model
% Compliance = Velocity / Pressure
% from displacement-stress eigenfunction of the model
% Input:
% frelist frequency array
% clist   acoustic wave / pressure drop horizontal speed, can be a scalar or a vector in the same size of frelist
% Vp, Vs, den P-wave velocity, S-wave velocity, density of each layer, N;
% z       the thickness of each layer except the halfspace, N-1;
%%%%
% This code is based on Aki&Richard (2002) Chapter 7.2
% Zongbo Xu 2022
%%%%

%% Check if the speed is of the same size as the frequecy list
if length(clist) == 1
    clist = ones(size(frelist))*clist;
end

pvcomp = ones(size(frelist));
pvcomph = ones(size(frelist));
tilt = zeros(size(frelist));
%% Loop through frequecies
for ifre = 1 : length(frelist)
    fre = frelist(ifre);
    k = fre * 2 * pi/clist(ifre);
%% eigenfunction values at the lower halfspace
omega = 2 * pi * fre;
y = sqrt(k^2 - omega^2/(Vp(end)^2));
v = sqrt(k^2 - omega^2/(Vs(end)^2));
mu = Vs(end)^2 * den(end);
% Equation 7.54
F = [Vp(end)*k Vs(end)*v; Vp(end)*y Vs(end)*k;-2*Vp(end)*mu*k*y -Vs(end)*mu*(k^2+v^2);-Vp(end)*mu*(k^2+v^2) -2*Vs(end)*mu*k*v];
Fwp = F * [exp(-y*sum(z)); 0]./omega;
Fws = F * [0; exp(-v*sum(z))]./omega; 
%% From the layer over the halfspace to the surface
fp = Fwp;
fs = Fws;
for ilayer = length(z):-1:1
P = pmatrix(Vp(ilayer),Vs(ilayer),den(ilayer),k,fre,-z(ilayer));
% Solution for P
fp = P * fp;
% for S
fs = P * fs;
end
%% Determine the ratio between P and S in the halfspace
% Zero shear stress
ratio = -fp(3)/fs(3);
f = fp  + fs* ratio;
%% Convert from displacement to velocity
% The very front -1 is due to different sign convesion between atmospheric
%    pressure and ground normal stress
% 2* pi *1i due to the ground displacement to velocity 
%    the IFFT convension in MATLAB is exp(i\omega t)
% an extra 1i in pvcomph is due to the pi/2 phase shift between Z and R 
%    components 
pvcomp(ifre) = -f(2)/f(4)* fre * 2*pi * (1i);
pvcomph(ifre) = -f(1)/f(4)* fre * 2*pi *(1i) * (1i);
% tilt
g = 3.72;
tilt(ifre) =  pvcomp(ifre)./ (2*pi*fre) ./ clist(ifre) * g;
end

end

function P=pmatrix(vp,vs,den,k,fre,dz)
% propagation matrix : P 
% Aki&Richard (2002) Equation 7.45
omega = 2 * pi * fre;
y = sqrt(k^2 - omega^2/(vp^2));
v = sqrt(k^2 - omega^2/(vs^2));
mu = vs^2 * den;
sinhy2 = sinh(y*dz/2);
sinhy = sinh(y*dz);
sinhv2 = sinh(v*dz/2);
sinhv = sinh(v*dz);
omden = omega^2*den;
P = zeros(4,4);
P(1,1) = 1 + 2 * mu /omden *( 2*k^2* sinhy2^2 - (k^2+v^2) * sinhv2^2 );
P(3,3) = P(1,1);
P(1,2) = k*mu/omden * ( (k^2 + v^2) * sinhy/y - 2*v*sinhv );
P(4,3) = - P(1,2);
P(1,3) = 1/omden * ( k^2*sinhy/y - v*sinhv );
P(1,4) = 2*k/omden * (sinhy2^2 - sinhv2^2 );
P(2,3) = -P(1,4);
P(2,1) = k*mu/omden * ( (k^2+v^2)*sinhv/v - 2*y*sinhy );
P(3,4) = -P(2,1);
P(2,2) = 1+2*mu/omden*( 2*k^2*sinhv2^2 - (k^2+v^2)*sinhy2^2 );
P(4,4) = P(2,2);
P(2,4) = 1/omden * ( k^2*sinhv/v - y*sinhy );
P(3,1) = mu^2/omden * (4*k^2*y*sinhy - (k^2+v^2)^2*sinhv/v );
P(3,2) = 2*mu^2*(k^2+v^2)*P(1,4);
P(4,1) = -P(3,2);
P(4,2) = mu^2/omden * ( 4*k^2*v*sinhv - (k^2+v^2)^2 * sinhy/y );
end
