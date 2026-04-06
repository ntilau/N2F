%%% Near field to far field transformation from a bounding box
clear; close all; clc;
addpath('../../matlabLib/');
arrayPos = buildArray(1, 9, .5, 1, .5);

points = 5;
smplRes = zeros(1,points);
nbrSmpls = smplRes;
error1 = smplRes;
error2 = smplRes;
errorM = smplRes;
for i=1:points
  smplRes(i) = 2^(-i); 
  [xMin, xMax, yMin, yMax, zMin, zMax, xPts, yPts, zPts] = ...
    getBoxDim(1, arrayPos, .5, smplRes(i), .5 );
  [boxPos, boxN, dS, mSize] = ...
    buildBox([1 1 1 1 1 1], xMin, xMax, yMin, yMax, zMin, zMax,...
    xPts, yPts, zPts, 1, 0, 0);
  [Rmag, NdotRV] = getBoxVectors(arrayPos, boxPos, boxN);
  excitPhasor = sf_excitations(1, arrayPos, 0, 0 );
  [psi, delPsi] = sf_nfSolver(1, excitPhasor, Rmag, NdotRV);
  dthetaFF = 1;
  thetaFF = deg2rad(-90:dthetaFF:90);
  phiFF = deg2rad([0 90]);
  fPsi = sf_nf2ffSolver(1, thetaFF, phiFF, boxPos, boxN, dS, psi, delPsi);
  fPsiRef = sf_directffSolver(1, thetaFF, phiFF, excitPhasor, arrayPos);

  nbrSmpls(i) = length(dS);
  error1(i) = getL1error(fPsi, fPsiRef);
  error2(i) = getL2error(fPsi, fPsiRef);
  errorM(i) = getMaxError(fPsi, fPsiRef);
end
%%
figProp= getFigureProperties();

figure;
% loglog(smplRes(i), error2,'*r','MarkerSize', figProp.ms);
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