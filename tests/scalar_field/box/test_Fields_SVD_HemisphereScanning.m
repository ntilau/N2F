%%% test of the nf2ff transformation from a sphere 
clear; clc; close all;
addpath('../../matlabLib');
tStart = tic;

arrayPos = buildArray(1, 3, .5, 5, .5);
[xMin, xMax, yMin, yMax, zMin, zMax, xPts, yPts, zPts] = ...
  getBoxDim(1, arrayPos, .5, .1, .5 );
[boxPos, boxN, dS, mSize] = ...
  buildBox([1 1 1 1 1 1], xMin, xMax, yMin, yMax, zMin, zMax,...
  xPts, yPts, zPts, 1, 0, 0);
[Rmag, NdotRV] = getBoxVectors(arrayPos, boxPos, boxN);

% ----- loop of testing
nbrMax=19;
trunc = 0;
colspanIdx = 0;
for i=1:nbrMax
  [spanningAnglesT, spanningAnglesP, testAnglesT, testAnglesP] = ...
  getSpiralingHelicoidalTrajectory(i,2*i, false, 1);
   
  spanPsi = zeros(size(boxPos,2),nbrMax);
  spanDelPsi = spanPsi;
  
  for j=1:length(spanningAnglesT)
    fprintf('--> angle : (%d)%.4g\n',j,angle(j));
    excitPhasor = sf_excitations(1, arrayPos, spanningAnglesT(j), spanningAnglesP(j) );
    [spanPsi(:,j), spanDelPsi(:,j)] = ...
      sf_nfSolver(1, excitPhasor, Rmag, NdotRV);
  end
  
  [uPsi, sPsi, vPsi] = svd(spanPsi,0);
  [uDelPsi, sDelPsi, vDelPsi] = svd(spanDelPsi,0);
  
  % --- SVD truncation
  if(trunc == 0)
    uPsiRed = uPsi;
    uDelPsiRed = uDelPsi;
  else
    uPsiRed = uPsi(:,1:trunc);
    uDelPsiRed = uDelPsi(:,1:trunc);
  end
  
  nError = zeros(1,length(testAnglesT));
  fError = nError;
  fRefError = nError;
  % --- Testing
  for p=1:length(testAnglesT)
    steering_t = testAnglesT(p);
    steering_p = testAnglesP(p);
    excitPhasor = sf_excitations(1, arrayPos, steering_t, steering_p );
    [tPsi, tDelPsi] = sf_nfSolver(1, excitPhasor, Rmag, NdotRV);

    thetaFF = deg2rad(-90:1:90);
    phiFF = deg2rad([0 90]);
    ftPsi = sf_nf2ffSolver(1, thetaFF, phiFF, boxPos, boxN, dS, ...
      tPsi, tDelPsi);
    ftPsiRef = sf_directffSolver(1, thetaFF, phiFF, ...
      excitPhasor, arrayPos);
    ftGainRef = sf_computeGain(ftPsiRef);
    tPsi = tPsi.';
    tDelPsi = tDelPsi.';

    aPsi = (uPsiRed* ((uPsiRed)' * tPsi));
    aDelPsi = (uDelPsiRed* ((uDelPsiRed)' * tDelPsi));

    faPsi = sf_nf2ffSolver(1, thetaFF, phiFF, boxPos, boxN, dS,...
      aPsi.', aDelPsi.');

    % --- Errors
    nError(p) = getL2error(aPsi+aDelPsi, tPsi+tDelPsi);
    fError(p) = getL2error(faPsi, ftPsi);
    fRefError(p) = getL2error(faPsi, ftPsiRef);
  end

  fErrorTot(i) = mean(fError); 

end
%%
figProp = getFigureProperties();
figure; semilogy(1:nbrMax,fErrorTot,'-r*', ...
  'LineWidth', figProp.lw, 'MarkerSize', figProp.ms);
v = axis;
axis([1 nbrMax v(3) v(4)]);
xlabel('Tested direction number', 'FontSize', figProp.fs);
ylabel('Relative error', 'FontSize', figProp.fs);
%printEPS('','errorPlanar');

fprintf('\nTotal computation time = %2.4g s\n', toc(tStart));
