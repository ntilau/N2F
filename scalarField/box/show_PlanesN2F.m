%%% Near field to far field transformation from the top and bottom planes 
%%% of a box encompassing a planar array of point sources
clear; close all; clc;
addpath('../../matlabLib/');

%%--- params
% c0 = 299792458;
% freq = 1e9;
lambda = 1; % c0/freq;
nbrElems_x = 3; % number of point sources on x direction
WLspacing_x = .5; % spacing between pointsources in wavelengths (wl.) x dir
nbrElems_y = 5; % analog on y
WLspacing_y = .5; % analog on y

%%--- nf Box params
WLranging = .5;
WLspacing = .1;
ext = 2.5;

%%--- Building Structure
arrayPos = buildArray(lambda, nbrElems_x, WLspacing_x, ...
  nbrElems_y, WLspacing_y);
[xMin, xMax, yMin, yMax, zMin, zMax, xPts, yPts, zPts] = ...
  getBoxDim(lambda, arrayPos, WLranging, WLspacing, ext );
[boxPos, boxN, dS, mSize] = ...
  buildBox([1 1 0 0 0 0], xMin, xMax, yMin, yMax, zMin, zMax,...
  xPts, yPts, zPts, 1, 0, 0);
% sf_plotArrayGeom(arrayPos);
[Rmag, NdotRV] = getBoxVectors(arrayPos, boxPos, boxN);

%%---
steering_t = 0; % [°]
steering_p = 0; % [°]
excitPhasor = sf_excitations(lambda, arrayPos, steering_t, steering_p );
[psi, delPsi] = sf_nfSolver(lambda, excitPhasor, Rmag, NdotRV);

% sf_plotBoxNF(mSize, boxPos, psi, '')
% printEPS('','boxfields')

%%--- nf2ff
dthetaFF = 1;
thetaFF = deg2rad(-90:dthetaFF:90);
phiFF = deg2rad([0 90]);
fPsi = sf_nf2ffSolver(lambda, thetaFF, phiFF, boxPos, boxN, dS, ...
  psi, delPsi);
fPsiRef = sf_directffSolver(lambda, thetaFF, phiFF, ...
  excitPhasor, arrayPos);
gain = sf_computeGain(fPsi);
gainRef = sf_computeGain(fPsiRef);
sf_plotFFCutPlanes(thetaFF, gain, thetaFF, gainRef, [1 2]);