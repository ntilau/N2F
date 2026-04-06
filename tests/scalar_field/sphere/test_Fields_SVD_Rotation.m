%%% low-rank approximation of the near fields with fields rotation around y
%%% axis
clear; clc; close all;
addpath('../../matlabLib');
tStart = tic();

dirToTest = 23.3251;
maximumSearch = true;
nbrElems_x=17;
arrayPos = buildArray(1, nbrElems_x, .5, 1, .5);
if maximumSearch
  ext = 1;
else
  ext = .5;
end
radius = getSphRadius(1, arrayPos, ext);
[spherePos, dS, thetaNF, phiNF, mSize] = buildSphere(radius, .1, ...
  3, 3, 1, dirToTest, [0 1 0]);
[Rmag, NdotRV, n] = getSphVectors(arrayPos, spherePos);
excitPhasor = sf_excitations(1, arrayPos, dirToTest, 0);
[tPsi, tDelPsi] = sf_nfSolver(1, excitPhasor, Rmag, NdotRV);
%% ----- reference for the tested angle
thetaFF = deg2rad(-90:1:90);
phiFF = deg2rad([0 90]);
ftPsi = sf_nf2ffSolver(1, thetaFF, phiFF, spherePos, n, dS, ...
  tPsi, tDelPsi);
ftPsiRef = sf_directffSolver(1, thetaFF, phiFF, ...
  excitPhasor, arrayPos);
%% ----- SVD testing
tPsi = tPsi.';
tDelPsi = tDelPsi.';
% sf_plotSphNFSolid(mSize, thetaNF, phiNF, tPsi, '');

range = 90;
maxNbrVects = 33;
[angles, nbrVectors] = getSpanningAngles(maxNbrVects, range);
nbrVectors = 1:nbrVectors(length(nbrVectors));

% for stAng=1:16
% dirToTest = angles(stAng) + (angles(stAng+1)-angles(stAng))/2;

spanPsi = zeros(size(spherePos,2),length(angles));
spanDelPsi = spanPsi;
nError = zeros(size(angles));
fError = nError;
fRefError = nError;

% ----- loop of testing

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

    [spherePos, dS, thetaNF, phiNF, mSize] = buildSphere(radius, .1, ...
      3, 3, 1, angle(j), [0 1 0]);
    [Rmag, NdotRV] = getSphVectors(arrayPos, spherePos);
    excitPhasor = sf_excitations(1, arrayPos, angle(j), 0 );
    [spanPsi(:,j), spanDelPsi(:,j)] = ...
      sf_nfSolver(1, excitPhasor, Rmag, NdotRV);
    
    if maximumSearch
      yIdx = find((abs(spherePos(2,:)) < .1) & (spherePos(3,:)>=0));
      psiTmp = spanPsi(:,j);
      [val, idx] = max(abs(psiTmp(yIdx)));
      idx = yIdx(idx);

      [p,t,r] = cart2sph(spherePos(1,idx), spherePos(2,idx), spherePos(3,idx));
      t=(pi/2-t)*180/pi;
%       if t>90
%         t=t-90;
%       end

      [spherePos, dS, thetaNF, phiNF, mSize] = buildSphere(radius, .1, ...
        3, 3, 1, t, [0 1 0]);
      [Rmag, NdotRV] = getSphVectors(arrayPos, spherePos);

      [spanPsi(:,j), spanDelPsi(:,j)] = ...
        sf_nfSolver(1, excitPhasor, Rmag, NdotRV);
    end
    sf_plotSphNFSolid(mSize, thetaNF, phiNF, spanPsi(:,j), '');
    pause(1)
  end
  close all;
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
  
  [spherePos, dS, thetaNF, phiNF, mSize] = buildSphere(radius, .1, ...
    3, 3, 1, dirToTest, [0 1 0]);
  [Rmag, NdotRV, n] = getSphVectors(arrayPos, spherePos);

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
% printEPS('',['SVDerrorExpRotated']);

fprintf('\nTotal computation time = %2.4g s\n', toc(tStart));
