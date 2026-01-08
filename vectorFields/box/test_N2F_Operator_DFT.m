%%% check for Vector N2F implementation with DFT-truncation of the
%%% operators

clear; clc; close all;
addpath('../../matlabLib');
tStart = tic;

% ----- load electrical parameters
[z0,k0,lambda0] = getFreeSpaceElectricalParams(1e9 , 0);

% ----- build Array
nbrElems_x = 3;
WLspacing_x = .5;
nbrElems_y = 5;
WLspacing_y = .5;
tJ=90; pJ=270; tM=90; pM=0; % [°] spherical components of Huygens' sources
Jmag=1; Mmag=z0; % [A/m] [V/m]
arrayPos = buildArray(lambda0, nbrElems_x, WLspacing_x, ...
  nbrElems_y, WLspacing_y);
steering_t = 0; % [°]
steering_p = 0; % [°]
[J,M] = vf_excitations(k0, arrayPos, Jmag, tJ, pJ,...
  Mmag, tM, pM, steering_t, steering_p);

% ----- build Box
WLranging = .5;
WLspacing = .1;
ext = .5;
[xMin, xMax, yMin, yMax, zMin, zMax, xPts, yPts, zPts] = ...
  getBoxDim(lambda0, arrayPos, WLranging, WLspacing, ext );
[boxPos, boxN, dS, mSize] = ...
  buildBox([1 1 1 1 1 1], xMin, xMax, yMin, yMax, zMin, zMax,...
  xPts, yPts, zPts, 1, 0, 0);

% ----- near Fields
[E, H, S, Pr] = vf_nfSolver(k0, z0, arrayPos, boxPos, boxN, dS, J, M);
nf = [H(:);E(:)];

% ----- far Fields
sizetref = 200;
fact = 2;
sizet = fact*sizetref;
nbrCoeffs = 21; % must be odd (if not the superio )
[phiFF, thetaFF] = meshgrid(deg2rad([0 90]), ...
  linspace(0, 2*pi*(sizetref-1)/sizetref, sizetref));
[Opt,Opp] =  vf_n2fOpFieldsFFT(k0, z0, boxPos, boxN, dS,...
  thetaFF, phiFF, nbrCoeffs);
[ifftOp, phiFFT, thetaFFT] = getLookAnglePoly(phiFF(1,:), sizet, ...
  sizetref, nbrCoeffs);
EtFF = zeros(sizet, size(phiFF,2));
EpFF = zeros(sizet, size(phiFF,2));
for i=1:size(phiFF,2)
  EtFF(:,i) = ifftOp*(Opt(:,:,i)*nf);
  EpFF(:,i) = ifftOp*(Opp(:,:,i)*nf);
end
[EtFFref, EpFFref] = vf_directffSolver(k0, z0, arrayPos, ...
  J, M, thetaFF, phiFF);

% ----- compute gain
[gaint, gainp] = vf_computeGain(z0, EtFF, EpFF, Pr);
[gaintRef, gainpRef] = vf_computeGain(z0, EtFFref, EpFFref, Pr);

% ----- compute error
errorN2F = getL2error(EtFF(1:fact:end,:)+EpFF(1:fact:end,:), EtFFref+EpFFref);

% ----- plot pattern
vf_plotFFPolarCutPlanes(figure(), 60, 10, ...
  thetaFFT, gaint, gainp,'N2F', thetaFF(1:3:end,:), ...
  gaintRef(1:3:end,:), gainpRef(1:3:end,:),'Direct', '\bf N2F Pattern');

fprintf('#> Total computation time : %g s.\n',toc(tStart))