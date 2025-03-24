function Atten = MarsAtmosphereAttenuation_BC(fre,T)
%% Acoustic wave attenuation in Martian atmosphere
% Output: 
%   Atten, the attnuation (column 2) varies with frequency (column 1)
% based on Bass and Chambers (2001)

if nargin == 1
    T = 200; %K
end

P = 670; % Pa
%T = 200; % K

if size(fre,1) == 1
    fre = fre.';
end

R = 8.3144621E3; % J/kmol/K
% S en K
S = 217;
% beta en kg/(m*s*(K**(1/2)))
beta = 1.49e-6;
% RM=R/44.0
theta = 960;
%Cvtrans = 3/2*R;
%Cvrot = 0.98*R;
% Cvrot=R
% this is for only CO2
Cv0 = (5/2)*R;
Cp0 = (7/2)*R;
gamma0 = 7/5;
Cvinf = (5/2)*R;
Cpinf = (7/2)*R;
% Mars air content in percentage
Xco2  =  0.953;
Xn2   =  0.027;
Xo2 = 0.0013;      % Changed from bass and chambers from 0.02
Xarg = 0.016;    % Added from Bass and Chambers
Xh2o = 0;%1 - (Xco2+Xn2+Xo2+Xarg);
M = 43.5; % kg/kmol


% CO2
MU = beta*sqrt(T)/(1+S/T); % kg/m
Cj = 2 * R * Xco2 * (theta/T)^2*exp(-theta/T)/(1-exp(-theta/T))^2;
Cv = R*3/2 + 0.98*R + Cj;
%%
%% Update from Raphael
% C^\prime_j
%Cprime = R*(2*Xco2.*((theta./T).^2)).*exp(-theta./T)./((1-exp(-theta./T)).^2);
% Cpinf = Cp0 - Cj;
% Cvinf = Cv0 - Cj;
%%
% Mayer's relation 
Cp = Cv + R; % only at 200K, not the case for 300K
gamma  = Cp/Cv;
c = sqrt(gamma*R*T/M);
kappa = (15*R*MU/4)*(4*Cv/(15*R)+3/5); % J*kg/kmol/K/m/s
%% Classic absorption
% viscosity and thermal conduction
Ac=((2*(pi*fre).^2)/(gamma*P*c)) * ((4/3)*MU+(gamma-1)*kappa/(gamma*Cv));
%% Rotation
Zrot=61.1*exp(-16.8/(T^(1/3)));
Ar=((2*(pi*fre).^2)/(gamma0*P*c)) * MU*gamma0*(gamma0-1)*R/(1.25*Cp0)*Zrot;
%% Vibration
kco2 = 0.219*P/MU*exp(-60.75/T^(1/3));
kn2  = 1.44*P/MU *exp(-78.29/T^(1/3));
karg = kn2;
kh2o = 6e-2*P/MU;
kk   = Xco2*kco2 +Xn2*kn2+ Xarg*karg +Xh2o*kh2o;
%
tauvt = 1/(kk *(1-exp(-theta/T)));
tauvs = (Cpinf/Cp) *tauvt;
fr =1/(2*pi*tauvs);
ss =  Cj*R/(Cpinf*(Cvinf+Cj));
Av = pi*ss/c * (fre.^2/fr)./(1+(fre/fr).^2);

Atten = zeros(length(fre),2);
Atten(:,1) = fre;
Atten(:,2) = Ac + Ar + Av;
end
