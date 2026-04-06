%%% test of the nf2ff transformation from a sphere 
clear; clc; close all;
addpath('../../matlabLib');

arrayPos = buildArray(1, 5, .5, 1, .5);
radius = getSphRadius(1, arrayPos, .5);
[spherePos, dS, thetaNF, phiNF, mSize] = buildSphere(radius, .1, 3, 3, 1);
[Rmag, NdotRV, n] = getSphVectors(arrayPos, spherePos);
excitPhasor = sf_excitations(1, arrayPos, 0, 0);
[psi, delPsi] = sf_nfSolver(1, excitPhasor, Rmag, NdotRV);
Psi = psi.';
DelPsi = delPsi.';
%% ff
thetaFF = deg2rad(-90:2:90);
phiFF = 0;
[A, B] = sf_nf2ffOperator(1, thetaFF, phiFF, ...
  spherePos, n, dS);
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
printEPS('','operatorSVD');