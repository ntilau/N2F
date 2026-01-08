%%% low-rank approximation - spanning space angles selection for low
%%% average error
clear; close all; clc;
addpath('../../matlabLib/');
tStart = tic;

arrayPos = buildArray(1, 21, .5, 5, .5);
[xMin, xMax, yMin, yMax, zMin, zMax, xPts, yPts, zPts] = ...
  getBoxDim(1, arrayPos, .5, .1, .5 );
[boxPos, boxN, dS, mSize] = ...
  buildBox([1 1 1 1 1 1], xMin, xMax, yMin, yMax, zMin, zMax,...
  xPts, yPts, zPts, 1, 0, 0);
[Rmag, NdotRV] = getBoxVectors(arrayPos, boxPos, boxN);

%%-- scan directions are selected closer to the broadside direction with an
%%-- exponential law
nbrVectors = 17;
coeff = 1.1;
angles = 90/(coeff^(nbrVectors-1)-1)* (coeff.^((1:nbrVectors)-1)-1);
trunc = nbrVectors;
%%-
nErrorStr = zeros(1,(nbrVectors-1));
fErrorStr = nErrorStr;
for stAng=1:(nbrVectors-1)
  
  % select test angle in the middle of the chosen spanning angles
  dirToTest = angles(stAng) + (angles(stAng+1)-angles(stAng))/2;
  excitPhasor = sf_excitations(1, arrayPos, dirToTest, 0 );
  [tPsi, tDelPsi] = sf_nfSolver(1, excitPhasor, Rmag, NdotRV);
  
  thetaFF = deg2rad(-90:.5:90);
  phiFF = deg2rad([0 90]);
  ftPsi = sf_nf2ffSolver(1, thetaFF, phiFF, boxPos, boxN, dS, ...
    tPsi, tDelPsi);
  ftPsiRef = sf_directffSolver(1, thetaFF, phiFF, ...
    excitPhasor, arrayPos);
  ftGainRef = sf_computeGain(ftPsiRef);

  spanPsi = zeros(size(boxPos,2),length(angles));
  spanDelPsi = spanPsi;
  nError = zeros(size(angles));
  fError = nError;
  fRefError = nError;

  % ----- loop of testing
  tPsi = tPsi.';
  tDelPsi = tDelPsi.';

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
      excitPhasor = sf_excitations(1, arrayPos, angle(j), 0 );
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

    faPsi = sf_nf2ffSolver(1, thetaFF, phiFF, boxPos, boxN, dS,...
      aPsi.', aDelPsi.');

    % --- Errors
    nError(colspanIdx) = getL2error(aPsi+aDelPsi, tPsi+tDelPsi);
    fError(colspanIdx) = getL2error(faPsi, ftPsi);
    fRefError(colspanIdx) = getL2error(faPsi, ftPsiRef);

    % --- Plot patterns
    if colspanIdx == 17 && false
      faGain = sf_computeGain(faPsi);
      ftGain = sf_computeGain(ftPsi);
      infos(1).legend1 = 'Tested N2F fields';
      infos(1).legend2 = 'Direct fields';
      infos(1).markers{1} = '-b';
      infos(1).markers{2} = '--r';
      sf_plotFFCutPlanes(thetaFF, faGain, thetaFF(:,1:3:end), ...
        ftGainRef(:,1:3:end), 1, infos);
%       printEPS('',['pattern_trunc',num2str(trunc),'_',num2str(colspanIdx)]);
    end

  end

  nErrorStr(stAng) = getL2error(aPsi+aDelPsi, tPsi+tDelPsi);
  fErrorStr(stAng) = getL2error(faPsi, ftPsi);
  
end
%%
figProp = getFigureProperties();
figure; semilogy(1:(nbrVectors-1),nErrorStr,'*', ...
  1:(nbrVectors-1),fErrorStr,'+r', ...
  'LineWidth', figProp.lw, ...
  'MarkerSize', figProp.ms);
v = axis;
axis([0 nbrVectors v(3) v(4)]);
xlabel('Tested direction number');
ylabel('Relative error');
legend('n.f. error', 'f.f. error');
%printEPS('','errorAvg');

fprintf('\nTotal computation time = %2.4g s\n', toc(tStart));
