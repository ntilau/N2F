%%% low-rank approximation of near fields: exponential angle selection
clear; clc; close all;
addpath('../../matlabLib');
tStart = tic;

dirToTest = 23.3251;
nbrElems_x = 9;
arrayPos = buildArray(1, nbrElems_x, .5, 5, .5);
radius = getSphRadius(1, arrayPos, .5);
[spherePos, dS, thetaNF, phiNF, mSize] = buildSphere(radius, .1, 3, 3, 1);
[Rmag, NdotRV, n] = getSphVectors(arrayPos, spherePos);
excitPhasor = sf_Excitations(1, arrayPos, dirToTest, 0);
[tPsi, tDelPsi] = sf_nfSolver(1, excitPhasor, Rmag, NdotRV);
%% ----- reference for the tested angle
thetaFF = deg2rad(-90:1:90);
phiFF = deg2rad(0);
ftPsi = sf_nf2ffSolver(1, thetaFF, phiFF, spherePos, n, dS, ...
  tPsi, tDelPsi);
ftPsiRef = sf_directffSolver(1, thetaFF, phiFF, ...
  excitPhasor, arrayPos);
%% ----- SVD testing
tPsi = tPsi.';
tDelPsi = tDelPsi.';

range = 90;
maxNbrVects = 33;
[angles, nbrVectors] = getSpanningAngles(maxNbrVects, range);
nbrVectors = 1:nbrVectors(length(nbrVectors));

spanPsi = zeros(size(spherePos,2),length(angles));
spanDelPsi = spanPsi;
nError = zeros(size(angles));
fError = nError;
fRefError = nError;
condNbr = zeros(size(nbrVectors));

% ----- loop of testing
flagSVD = false;
trunc = 0;
colspanIdx = 0;
initial = find(nbrVectors>=trunc, 1);
for i=initial:length(nbrVectors)
  angle = angles(1:nbrVectors(i));
  colspanIdx = colspanIdx + 1; 
  fprintf('++++++++++++++++++++++++++++++++++++++\n');
  fprintf('Nbr of Vectors = %g, Direction to test = %.8g\n',...
      nbrVectors(i), dirToTest);
  
  for j=colspanIdx:length(angle)
    fprintf('--> angle : (%d)%.4g\n',j,angle(j));
    excitPhasor = sf_Excitations(1, arrayPos, angle(j), 0 );
    [spanPsi(:,j), spanDelPsi(:,j)] = ...
      sf_nfSolver(1, excitPhasor, Rmag, NdotRV);
  end
  
  % --- SVD
  colspanIdx = length(angle);
  [uPsi, sPsi, vPsi] = svd(spanPsi(:,1:colspanIdx),0);
  [uDelPsi, sDelPsi, vDelPsi] = svd(spanDelPsi(:,1:colspanIdx),0);
  
  % --- SVD truncation
  if(trunc == 0)
    uPsiRed = uPsi;
    uDelPsiRed = uDelPsi;
  else
    uPsiRed = uPsi(:,1:trunc);
    uDelPsiRed = uDelPsi(:,1:trunc);
  end
  
  % --- Testing
  aPsi = (uPsiRed* ((uPsiRed)' * tPsi));
  aDelPsi = (uDelPsiRed* ((uDelPsiRed)' * tDelPsi));

  faPsi = sf_nf2ffSolver(1, thetaFF, phiFF, spherePos, n, dS,...
    aPsi.', aDelPsi.');
  
  % --- Errors
  nError(colspanIdx) = getL2error(aPsi+aDelPsi, tPsi+tDelPsi);
  fError(colspanIdx) = getL2error(faPsi, ftPsi);
  fRefError(colspanIdx) = getL2error(faPsi, ftPsiRef);

end


% --- Plot SVD errors
fprintf('Reference :\n');
refError = getL2error(ftPsi, ftPsiRef);
plotSVDerror(sPsi, nbrVectors, nError, fError, fRefError, nbrElems_x, false);
printEPS('','SVDerrorExp_33');

fprintf('\nTotal computation time = %2.4g s\n', toc(tStart));
