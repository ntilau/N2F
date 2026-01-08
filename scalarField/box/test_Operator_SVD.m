%%% test of the nf2ff transformation from a sphere 
clear all; clc; close all;
addpath('..\..\matlabLib');

arrayPos = buildArray(1, 5, .5, 1, .5);
[xMin, xMax, yMin, yMax, zMin, zMax, xPts, yPts, zPts] = ...
  getBoxDim(1, arrayPos, .5, .1, .5 );
[boxPos, boxN, dS, mSize] = ...
  buildBox([1 1 1 1 1 1], xMin, xMax, yMin, yMax, zMin, zMax,...
  xPts, yPts, zPts, 1, 0, 0);
[Rmag, NdotRV] = getBoxVectors(arrayPos, boxPos, boxN);
excitPhasor = sf_Excitations(1, arrayPos, 0, 0);
[psi, delPsi] = sf_nfSolver(1, excitPhasor, Rmag, NdotRV);
Psi = psi.';
DelPsi = delPsi.';
%% ff
thetaFF = deg2rad(-90:2:90);
phiFF = 0;
[A, B] = sf_nf2ffOperator(1, thetaFF, phiFF, ...
  boxPos, boxN, dS);
fPsiRef = sf_directffSolver(1, thetaFF, phiFF, ...
  excitPhasor, arrayPos);
Size = size(A);
%% plot singular values
[ua,sa,va] = svd(A.',0);
figProp = getFigureProperties();
semilogy(diag(sa),'o-b', 'LineWidth', figProp.lw, ...
  'MarkerSize', figProp.ms);
axis tight
xlabel('n', 'FontSize', figProp.fs);
ylabel('Singular values \sigma_n', 'FontSize', figProp.fs);
% printEPS('','operatorSVD');