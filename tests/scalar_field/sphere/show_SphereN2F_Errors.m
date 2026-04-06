%%% check the error committed for several sampling resolutions
clear; clc; close all;
addpath('../../matlabLib');

arrayPos = buildArray(1, 3, .5, 5, .5);

points = 5;
smplRes = zeros(1,points);
nbrSmpls = smplRes;
error1 = smplRes;
error2 = smplRes;
errorM = smplRes;
for i=1:points
  sphSmplRes = .5/2.^(i); 
  radius = getSphRadius(1, arrayPos, .5);
  [spherePos, dS, thetaNF, phiNF, matrixSize] = ...
    buildSphere(radius, sphSmplRes, 3, 3, 1);
  [Rmag, NdotRV, n] = getSphVectors(arrayPos, spherePos);
  excitPhasor = sf_excitations(1, arrayPos, 0, 0 );
  [psi, delPsi] = sf_nfSolver(1, excitPhasor, Rmag, NdotRV);
  dthetaFF = 1;
  thetaFF = deg2rad(-90:dthetaFF:90);
  phiFF = deg2rad([0 90]);
  fPsi = sf_nf2ffSolver(1, thetaFF, phiFF, spherePos, n, dS, psi, delPsi);
  fPsiRef = sf_directffSolver(1, thetaFF, phiFF, excitPhasor, arrayPos);

  smplRes(i)= sphSmplRes;
  nbrSmpls(i) = length(dS);
  error1(i) = getL1error(fPsi, fPsiRef);
  error2(i) = getL2error(fPsi, fPsiRef);
  errorM(i) = getMaxError(fPsi, fPsiRef);
end
%%
figProp= getFigureProperties();

figure;
% loglog(smplRes, error2,'*r','MarkerSize', figProp.ms);
loglog(nbrSmpls, error2,'*r','MarkerSize', figProp.ms);
xlabel('Number of samples', 'FontSize', figProp.fs);
ylabel('Relative error', 'FontSize', figProp.fs);
legend('\bf{L}_2','Location','SouthEast');

figure;
loglog(smplRes, error1,'+',smplRes, error2,'*r',smplRes, ...
  errorM,'ok','MarkerSize', figProp.ms);
xlabel('Sampling resolution [\lambda]', 'FontSize',figProp.fs);
ylabel('Relative error', 'FontSize',figProp.fs);
legend('\bf{L}_1','\bf{L}_2','\bf{L}_\infty','Location','SouthEast');


