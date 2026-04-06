%%% check for Vector N2F implementation

clear; clc; close all;
addpath('../../matlabLib');

% ----- load electrical parameters
[z0,k0,lambda0] = getFreeSpaceElectricalParams(1e9 , 0);

% ----- build Array
nbrElems_t = 3;
WLspacing_t = .5;
nbrElems_p = 5;
WLspacing_p = .5;
tJ=90; pJ=270; tM=90; pM=0; % [°] spherical components of Huygens' sources
Jmag=1; Mmag=z0; % [A/m] [V/m]
arrayPos = buildArray(lambda0, nbrElems_t, WLspacing_t, ...
  nbrElems_p, WLspacing_p);
steering_t = 0; % [°]
steering_p = 0; % [°]
[J,M] = vf_excitations(k0, arrayPos, Jmag, tJ, pJ,...
  Mmag, tM, pM, steering_t, steering_p);

% ----- build Box
WLranging = .5;
WLspacing = .05;
ext = .5;
[xMin, xMax, yMin, yMax, zMin, zMax, xPts, yPts, zPts] = ...
  getBoxDim(lambda0, arrayPos, WLranging, WLspacing, ext );
[boxPos, boxN, dS, mSize] = ...
  buildBox([1 1 1 1 1 1], xMin, xMax, yMin, yMax, zMin, zMax,...
  xPts, yPts, zPts, 1, 0, 0);

% ----- near Fields
[E, H, S, Pr] = vf_nfSolver(k0, z0, arrayPos, boxPos, boxN, dS, J, M);

% ----- far Fields
% [phiFF, thetaFF] = meshgrid((-180:2:180)*pi/180, (-180:2:180)*pi/180);
[phiFF, thetaFF] = meshgrid([0 90]*pi/180, (-180:2:180)*pi/180);
[EtFF, EpFF] =  vf_nf2ffSolver(k0, z0, boxPos, boxN, dS,...
  E, H, thetaFF, phiFF);
[EtFFref, EpFFref] =  vf_directffSolver(k0, z0, arrayPos, ...
  J, M, thetaFF, phiFF);

% ----- load FEKO.txt
FEKO = dlmread('FEKO.txt', ' ', 0, 0);
FEKOEtFF = [FEKO(1:181,4).*exp(1i*FEKO(1:181,5)*pi/180), ...
  FEKO(182:362,4).*exp(1i*FEKO(182:362,5)*pi/180)];
FEKOEpFF = [FEKO(1:181,6).*exp(1i*FEKO(1:181,7)*pi/180), ...
  FEKO(182:362,6).*exp(1i*FEKO(182:362,7)*pi/180)];

% ----- compute errors
errorFEKOref = getL2error(FEKOEtFF+FEKOEpFF, EtFFref+EpFFref);
errorFEKOn2f = getL2error(FEKOEtFF+FEKOEpFF, EtFF+EpFF);
errorN2F = getL2error(EtFF+EpFF, EtFFref+EpFFref);

% ----- compute gain
[gaint, gainp] = vf_computeGain(z0, EtFF, EpFF, Pr);
[gaintRef, gainpRef] = vf_computeGain(z0, EtFFref, EpFFref, Pr);
[gaintFEKO, gainpFEKO] = vf_computeGain(z0, FEKOEtFF, FEKOEpFF, Pr);

% ----- plot FEKO pattern
handle = figure();
vf_plotFFPolarCutPlanes(handle, 60, 10, ...
  thetaFF, gaintFEKO, gainpFEKO, 'FEKO',thetaFF(1:3:end,:), ...
  gaintRef(1:3:end,:), gainpRef(1:3:end,:),'Direct', '\bf FEKO pattern');

% ----- plot radiation solid
% vf_plotFF3d(figure(), thetaFF, phiFF, 60, gaint, gainp)
% printEPS('','radSolid');

% ----- plot pattern
% handle = figure();
handle = [figure(),figure()];
vf_plotFFPolarCutPlanes(handle, 60, 10, ...
  thetaFF, gaint, gainp,'N2F', thetaFF(1:3:end,:), ...
  gaintRef(1:3:end,:), gainpRef(1:3:end,:),'Direct', '\bf N2F Pattern');
% printEPS('','errorPlot');
% for i=1:size(handle,2)
%   pause(.5)
%   figure(handle(i));
%   pause(.1);
%   printEPS('',['test',num2str(i)]);
% end
