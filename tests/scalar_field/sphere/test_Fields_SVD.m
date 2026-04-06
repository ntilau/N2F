%%% low-rank approximation of the near fields
clear; clc; close all;
addpath('../../matlabLib');
tStart = tic;

dirToTest = 23.3251;
arrayPos = buildArray(1, 9, .5, 5, .5);
radius = getSphRadius(1, arrayPos, .5);
[spherePos, dS, thetaNF, phiNF, mSize] = buildSphere(radius, .1, 3, 3, 1);
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

range = 90;
maxNbrVects = 40;
[angles, nbrVectors] = getSpanningAngles(maxNbrVects, range);
% angles = rand(size(angles))*90;

% ---- loop of testing
trunc=0;
colspanIdx=1;
for i=1:length(nbrVectors)
  angle = angles(1:nbrVectors(i));
  fprintf('++++++++++++++++++++++++++++++++++++++\n');
  fprintf('Nbr of Vectors = %g, Direction to test = %.8g\n',...
      nbrVectors(i), dirToTest);
  fprintf('--> angles : ');

  for j=colspanIdx:length(angle)
    excitPhasor = sf_excitations(1, arrayPos, angle(j), 0 );
    [spanPsi(:,j), spanDelPsi(:,j)] = ...
      sf_nfSolver(1, excitPhasor, Rmag, NdotRV);
    fprintf('(%d)%.4g ',j,angle(j));
  end
  
  colspanIdx = length(angle) + 1;
  fprintf('\n');

  [uPsi, sPsi, vPsi] = svd(spanPsi,0);
  [uDelPsi, sDelPsi, vDelPsi] = svd(spanDelPsi,0);

  if(trunc == 0)
    uPsiRed = uPsi;
    uDelPsiRed = uDelPsi;
  else
    uPsiRed = uPsi(:,1:trunc);
    uDelPsiRed = uDelPsi(:,1:trunc);
  end

  aPsi = (uPsiRed* ((uPsiRed)' * tPsi));
  aDelPsi = (uDelPsiRed* ((uDelPsiRed)' * tDelPsi));
  
  faPsi = sf_nf2ffSolver(1, thetaFF, phiFF, spherePos, n, dS,...
    aPsi.', aDelPsi.');
  
  %% Errors computed at broadside
  nError(i) = getL2error(aPsi, tPsi);
  fError(i) = getL2error(faPsi, ftPsi);
end
%% Error plot
figure
subplot(1,2,1);
semilogy(nbrVectors, nError,'.-');
hold on
semilogy(nbrVectors, fError,'*-r');
axis tight;
legend('Rel NF error', 'Rel FF error','Location','NorthEast');
subplot(1,2,2);
semilogy(diag(sPsi),'.-');
axis tight;

fprintf('\nTotal computation time = %2.4g s\n', toc(tStart));
