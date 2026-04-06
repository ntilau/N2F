%%% Check of the error due to numerical integration of the near fields in the
%%% vector Huygen's principle
clear; clc; close all;
addpath('../../matlabLib');

%%--- electrical parameters
[z0,k0,lambda0] = getFreeSpaceElectricalParams(1e9 , 0);

%%--- array parameters
nbrElems_x = 3; % number of point sources on x direction
WLspacing_x = .5; % spacing between pointsources in wavelengths (wl.) x dir
nbrElems_y = 5; % analog on y
WLspacing_y = .5; % analog on y
tJ=90; pJ=270; tM=90; pM=0; % Huygens' sources
Jmag=1; Mmag=z0; % [A/m] [V/m]

%%--- nf Box params
for i=2.^(0:6)
  
  WLranging = .5;
  WLspacing = .5/i;
  fprintf('-----------------\nWLspacing = %2.4g\n-----------------\n', ...
    WLspacing);
  ext = .5;

  arrayPos = buildArray(lambda0, nbrElems_x, WLspacing_x, ...
    nbrElems_y, WLspacing_y);
  steering_x = 0; % [°]
  steering_y = 0; % [°]
  [J,M] = vf_excitations(k0, arrayPos, Jmag, tJ, pJ,...
    Mmag, tM, pM, steering_x, steering_y);
  [xMin, xMax, yMin, yMax, zMin, zMax, xPts, yPts, zPts] = ...
    getBoxDim(lambda0, arrayPos, WLranging, WLspacing, ext );
  [boxPos, boxN, dS, mSize] = ...
    buildBox([1 1 1 1 1 1], xMin, xMax, yMin, yMax, zMin, zMax,...
    xPts, yPts, zPts, 1, 0, 0);
  [E, H, S, Pr] = vf_nfSolver(k0, z0, arrayPos, boxPos, boxN, dS, J, M);
  %%--- ff parameters
  % [phiFF, thetaFF] = meshgrid((-180:2:180)*pi/180, (-180:2:180)*pi/180);
  [phiFF, thetaFF] = meshgrid([0 90]*pi/180, (-180:2:180)*pi/180);
  [EtFF, EpFF] =  vf_nf2ffSolver(k0, z0, boxPos, boxN, dS,...
    E, H, thetaFF, phiFF);
  [EtFFref, EpFFref] =  vf_directffSolver(k0, z0, arrayPos, ...
    J, M, thetaFF, phiFF);

  boxSampling(i)= WLspacing;
  boxerror2(i) = getL2error(EtFF+EpFF, EtFFref+EpFFref);
end
%%
figProp= getFigureProperties();
figure;
loglog(boxSampling, boxerror2,'*r','MarkerSize',figProp.ms);
xlabel('Sampling resolution [\lambda]', 'FontSize',figProp.fs);
ylabel('Relative error', 'FontSize',figProp.fs);
printEPS('', 'boxSamplingVect');