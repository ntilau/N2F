%%% low-rank approximation of near fields: exponential angle selection and
%%% truncation of the basis
clear all; clc; close all;
addpath('..\..\matlabLib');
tStart = tic;

dirToTest = 23.3251;
nbrElems_x = 17;
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
nbrVec = [17 33];
for nbrVecIdx = 1:length(nbrVec)
  maxNbrVects = nbrVec(nbrVecIdx);
  [angles, nbrVectors] = getSpanningangles(maxNbrVects, range);

  spanPsi = zeros(size(spherePos,2),length(angles));
  spanDelPsi = spanPsi;
  nError = zeros(size(angles));
  fError = nError;
  fRefError = nError;

  truncLevel = 10:17;
  nError_trunc = zeros(size(truncLevel));
  fError_trunc = nError_trunc;
  diagPsi = nError_trunc;

  % ----- loop of testing
  colspanIdx = 0;
  initial = length(nbrVectors);%find(nbrVectors>=trunc, 1);
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
    for tr=1:length(truncLevel)
      trunc = truncLevel(length(truncLevel)-tr+1);

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
      nError(colspanIdx) = getL2error(aPsi, tPsi);
      fError(colspanIdx) = getL2error(faPsi, ftPsi);
      fRefError(colspanIdx) = getL2error(faPsi, ftPsiRef);

      nError_trunc(tr) = norm(aPsi-tPsi);
      fError_trunc(tr) = norm(faPsi-ftPsi);

    end
  end
% write to file for table
  a(:,2*(nbrVecIdx-1)+(1:2)) = [nError_trunc.', fError_trunc.' ];
end

dlmwrite('truncationErrorSigma1.txt',a, 'precision', 4);

fprintf('\nTotal computation time = %2.4g s\n', toc(tStart));
