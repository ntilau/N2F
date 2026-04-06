%%% test of the nf2ff transformation from a triangulated patches sphere 
%%% surrounding a planar array
clear; clc; close all;
addpath('../../matlabLib');

arrayPos = buildArray(1, 3, .5, 5, .5);
%%--- sphere
ext = .5; % gap between array and sphere in wavelengths
triOrder = 3; % mesh refinement order
plotTriSphere = true; % plot the sphere
radius = getSphRadius(1, arrayPos, ext);
[spherePos, dS, n] = buildTriSphere(radius, triOrder, plotTriSphere);
[Rmag, NdotRV] = getSphVectors(arrayPos, spherePos);
excitPhasor = sf_excitations(1, arrayPos, 15, 0);
[psi, delPsi]=  sf_nfSolver(1, excitPhasor, Rmag, NdotRV);
%%--- nf2ff
dtheta = .5; % ff pattern resolution [°]
theta = deg2rad(-90:dtheta:90);
phi = deg2rad([0 90]);
fPsi = sf_nf2ffSolver(1, theta, phi, spherePos, n, dS, ...
  psi, delPsi);
fPsiRef = sf_directffSolver(1, theta, phi, ...
  excitPhasor, arrayPos);
%%--- pattern plot
gain = sf_computeGain(fPsi);
refGain = sf_computeGain(fPsiRef);
sf_plotFFCutPlanes(theta, gain, theta, refGain, [1 2]);
error = getL2error(fPsi,fPsiRef);