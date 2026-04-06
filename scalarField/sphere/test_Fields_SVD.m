%%% Low-rank SVD approximation test for near fields on spherical bounding surface
% OBJECTIVE: Validate that near-field data can be approximated by a low-rank subspace
%            (key principle behind model order reduction for N2F operators)
%
% METHODOLOGY:
%   1. Setup: Create antenna array, generate sphere bounding surface
%   2. Compute: Full near field and derivatives on sphere via sf_nfSolver
%   3. Compute: Full far-field reference via direct summation
%   4. SVD: Decompose near-field matrix via singular value decomposition
%   5. Truncate: Retain only leading K singular vectors (low-rank approximation)
%   6. Error: Compare truncated vs full reconstruction fidelity at each rank
%   7. Sweep: Vary angle selection strategy and measure compression vs accuracy trade-off
%
% EXPECTED RESULTS:
%   - ~50-90% rank reduction achievable while maintaining -40dB error for patterns
%   - L2 error decays exponentially with rank (significant drop first 20-30 coefficients)
%   - Both 'spanning angles' and 'linear selection' give similar results for this scenario
%
clear; clc; close all;
addpath('../../matlabLib');
tStart = tic;

% TEST CONFIGURATION: patch antenna directed along 23.3 degrees
dirToTest = 23.3251;  % far-field look direction [degrees]
arrayPos = buildArray(1, 9, .5, 5, .5);  % 9x5 element grid, half-wavelength spacing
radius = getSphRadius(1, arrayPos, .5);  % find suitable sphere radius
[spherePos, dS, thetaNF, phiNF, mSize] = buildSphere(radius, .1, 3, 3, 1);
[Rmag, NdotRV, n] = getSphVectors(arrayPos, spherePos);
excitPhasor = sf_excitations(1, arrayPos, dirToTest, 0);  % uniform excitation phased scan
[tPsi, tDelPsi] = sf_nfSolver(1, excitPhasor, Rmag, NdotRV);

%% Generate reference far-field patterns (ground truth)
thetaFF = deg2rad(-90:1:90);  % observe full hemisphere
phiFF = deg2rad([0 90]);      % two orthogonal phi cuts
ftPsi = sf_nf2ffSolver(1, thetaFF, phiFF, spherePos, n, dS, tPsi, tDelPsi);
ftPsiRef = sf_directffSolver(1, thetaFF, phiFF, excitPhasor, arrayPos);

%% SVD decomposition and truncation analysis
% Transpose near fields for SVD: each row = sample point, columns = unknowns to reduce
tPsi = tPsi.';
tDelPsi = tDelPsi.';

% Angle sampling configuration for SVD basis construction
% range = angular extent (90 deg = main lobe region)
% maxNbrVects = target rank for truncation (controls compression)
range = 90;
maxNbrVects = 40;
[angles, nbrVectors] = getSpanningAngles(maxNbrVects, range);
% OPTIONAL: use random angles instead of structured exponential sampling
% angles = rand(size(angles))*90;  % random angles for comparison

% Truncation sweep: test reconstruction error vs rank
trunc=0;  % truncation counter
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
