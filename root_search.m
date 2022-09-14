function [fphv,fgv] = root_search(alpha,w,rho,thickness,nroot)
% Search the root of Equation 12: \partial z f = A f
%% Input
% alpha: the acoustic-wave velocity in each layer, from the ground surface
%        to the top halfspace;
% w:     in-line wind speed in each layer;
% rho:   air density in each layer;
% thickness: thickness of each alyer;
% nroot: root number at each frequency. 
% nroot = 1, fundamental mode
% nroot = 2, fundamental and 1st-higher mode
%% Output
% fphv: frequency and the phase velocity
% fgv: frequency and the group velocity
%%
if nargin == 4
   nroot = 1;
end

%% velocity search parameters
dc = 9E-4;
middv = dc/15;%3E-4
fres = 0.1:0.1:3;
fphv = zeros(length(fres),1+nroot);
fphv(:,1) = fres';
fgv = zeros(length(fres),1+nroot);
fgv(:,1) = fres';
alphamin = min(alpha+w);
alphamax = max(alpha+w);
%% Loop through frequencies to find phase velocity
for ifre = 1 : length(fres) 
    fre = fres(ifre);
    omega = 2 * pi *fre;
    iroot = 1;
    for c = alphamin : dc : alphamax
        a = c; b = c + dc;
        funa = PropMatW(alpha,w,rho,thickness,omega,a);
        funb = PropMatW(alpha,w,rho,thickness,omega,b);
        if real(funa)*real(funb) > 0
            continue;
        else
            while abs(a-b) >= middv
              midab = (a+b)/2;
              funmid = PropMatW(alpha,w,rho,thickness,omega,midab);
              if funa * funmid >= 0
                  a = midab;
                  funa = funmid;
              else
                  b = midab;
              end
            end
            if iroot <= nroot
               fphv(ifre,1 + iroot) = (a+b)/2;
               iroot = iroot + 1;
            else
               break
            end
        end
    end
            
end

% phase velocity to group velocity
for iroot = 1 : nroot
    id = fphv(:,iroot+1)>0;
    domega = 2*pi*gradient(fphv(id,1));
    dk     = gradient(2*pi*fphv(id,1)./fphv(id,iroot+1));
    fgv(id,iroot+1) = domega ./ dk;
end

end
