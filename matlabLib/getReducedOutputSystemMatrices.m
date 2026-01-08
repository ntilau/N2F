function [romCt, romCp, romF, nbrSmpls, nbrCoeffs] = ... 
  getReducedOutputSystemMatrices( nbrCoeffs, ...
  nbrSmpls, k0, z0, boxPos, boxN, dS, phi, Q)

% nbrCoeffs = (nbrCoeffs-1)/2;
disp('#> Loading Output System matrices ...');
comp2Ext = mmread('lte_fileset\comp2Ext1.mm');
num2UnNum = mmread('lte_fileset\num2UnNum1.mm');
mFE = mmread('lte_fileset\functionalE.mm');
mFH = mmread('lte_fileset\functionalH.mm');
fieldFunctional  = [mFH;mFE];
% coeffLevels = 1;

fprintf('#> Computing the n2f Fields operators ... ')
tic();
% while minCoeffsLevel < coeffLevels
%   nbrSmpls = floor(nbrCoeffs * overSampling);
  coeffLevels = 0;
  for i=1:size(phi,2)
    [phiFF,thetaFF] = meshgrid(phi(1,i),...
      linspace(0,2*pi*(nbrSmpls-1)/nbrSmpls,nbrSmpls));
    [Opt, Opp] = n2fOpFieldsFFT(k0, z0, boxPos, boxN, dS, ...
        thetaFF, phiFF, nbrCoeffs);
    level = max(max(max(abs(Opt(nbrCoeffs+1,:,:)))), ...
      max(max(abs(Opp(nbrCoeffs+1,:,:)))));
    coeffLevels = max(coeffLevels,level);
  end
%   nbrCoeffs = nbrCoeffs + 1;
% end
fprintf('Elapsed %2.4g s.\n', toc());

% nbrCoeffs = nbrCoeffs - 1;
% nbrSmpls = floor(nbrCoeffs * overSampling);% floor(nbrCoeffs * overSampling);
fprintf('##> Coeff levels = %2.4g, nbrCoeffs = %d\n', ...
  coeffLevels, nbrCoeffs*2+1);

fprintf('#> Computing the Reduced Order n2f Fields operators ...\n')
tic();
romCt = zeros(nbrCoeffs*2+1,size(Q,2),size(phi,2));
romCp = romCt;
for i=1:size(phi,2)
  fprintf('##> Phi = %2.4g°\n', phi(1,i)*180/pi);
  [phiFF,thetaFF] = meshgrid(phi(1,i),...
    linspace(0,2*pi*(nbrSmpls-1)/nbrSmpls,nbrSmpls));
  [Opt, Opp] = n2fOpFieldsFFT(k0, z0, boxPos, boxN, dS, ...
      thetaFF, phiFF, nbrCoeffs);
  romCt(:,:,i) = Opt*fieldFunctional*num2UnNum*comp2Ext*Q;
  romCp(:,:,i) = Opp*fieldFunctional*num2UnNum*comp2Ext*Q;
  
end
romF = fieldFunctional*num2UnNum*comp2Ext*Q;
fprintf('Elapsed %2.4g s.\n', toc());
