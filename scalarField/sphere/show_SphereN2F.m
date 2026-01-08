%%% test of the near field to far field transformation from a bounding 
%%% sphere 
clear; clc; close all;
addpath('../../matlabLib');

%%--- params for planar array of point sources on XY plane
lambda = 1;
nbrElems_x = 5; % number of point sources on x direction
WLspacing_x = .5; % spacing between pointsources in wavelengths (wl.) x dir
nbrElems_y = 5; % analog on y
WLspacing_y = .5; % analog on y

%%--- nf Sphere params
flag = 0; % 0: dThetaNF-dPhiNF  1: sampling resolution 2: for nf plots 
          % n.b. 2: induces errors in the n2f computation
sphSmplRes = .05*lambda; % sampling resolution on the sphere -> derives the
                        % angle resolution dThetaNF and dPhiNF
dThetaNF = 5;
dPhiNF = 5;
ext = .5; % minimum distance in wl. between the sphere and the array

%%--- building structure
arrayPos = buildArray(lambda, nbrElems_x, WLspacing_x, ...
  nbrElems_y, WLspacing_y);
radius = getSphRadius(lambda, arrayPos, ext);
[spherePos, dS, thetaNF, phiNF, matrixSize] = ...
  buildSphere(radius, sphSmplRes, dThetaNF, dPhiNF, flag);
[Rmag, NdotRV, n] = getSphVectors(arrayPos, spherePos);
plotSphGeom(matrixSize, arrayPos, spherePos);

%%--- nf
steering_t=45;
steering_p=0;
excitPhasor = sf_excitations(lambda, arrayPos, steering_t, steering_p);
[psi, delPsi] = sf_nfSolver(lambda, excitPhasor, Rmag, NdotRV);
% sf_plotSphNF(matrixSize, thetaNF, phiNF, psi.', '', true);
% sf_plotSphNFSolid(matrixSize, thetaNF, phiNF, psi.', '');

%%--- nf2ff
dthetaFF = .5; % ff pattern resolution [°]
thetaFF = deg2rad(-90:dthetaFF:90);
phiFF = deg2rad([0 90]);
fPsi = sf_nf2ffSolver(lambda, thetaFF, phiFF, spherePos, n, dS, ...
  psi, delPsi);
fPsiRef = sf_directffSolver(lambda, thetaFF, phiFF, ...
  excitPhasor, arrayPos);

%%--- plots
gain = sf_computeGain(fPsi);
refGain = sf_computeGain(fPsiRef);
planes = [1 2];
infos(1).title = {' Pattern in \phi = 0°'};
infos(2).title = {' Pattern in \phi = 90°'};
infos(1).legend1 = 'N2F';
infos(1).legend2 = 'Direct';
infos(2).legend1 = 'N2F';
infos(2).legend2 = 'Direct';
% filename{1} = 'sph3x5_Phi0_15_0';
% filename{2} = 'sph3x5_Phi90_15_0';
sf_plotFFCutPlanes(thetaFF, gain, thetaFF, refGain, planes, infos);