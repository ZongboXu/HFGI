function AcousticWavGuideWind_MARS_FundFirst

HOME = './';

%thickness = [20,20,40,55,55,115,115,180,200,200,200,300];
%altitude = [20,40,80,135,190,305,420,600,800,1000,1200,1500];

nroot = 3;

% Table S3
model = load('S0986c_ATM.dat');
thickness = model(1:end-1,1);
alpha = model(:,2);
w = model(:,3);
rho = model(:,4);
% 
% fundamental + 1st higher mode
[fphv,fgv] = root_search(alpha,w,rho,thickness,nroot);

norm_P  = zeros(length(fphv),nroot+1);
norm_P(:,1) = fphv(:,1);
for iroot = 1 : nroot
    for ifre = 1 : length(fphv(:,1))
        if fphv(ifre,iroot+1) > 0
            fprintf('%i mode \n',iroot-1);
            omega = fphv(ifre,1) * 2 * pi;
            norm_P(ifre,1+iroot) = Intergral(alpha,w,rho,thickness,omega,fphv(ifre,iroot+1));
        end
    end
end

% save('DispPhv.dat','fphv','-ascii');
% save('DispGv.dat','fgv','-ascii');
% save('NormP.dat','norm_P','-ascii');

end