% Computes the sphere sampling angles for plot purposes. As the initial and
% final angles coincides, this introduces an error in the near field to far
% field transformation
%
% [theta, phi, matrixSize] = getSphSmplAngles(dTheta, dPhi)
%
% IN: dTheta, dPhi = theta and phi resolutions in degrees
% 
% OUT: theta, phi = column vectors of the angles
%      matrixSize = size of the angles meshgrid for plots
%
% Laurent Ntibarikure
function [theta, phi, matrixSize] = getSphSmplAnglesForPlots(dTheta, dPhi)

[theta, phi] = meshgrid(deg2rad(0:dTheta:180),deg2rad(0:dPhi:360));
matrixSize = size(theta);
theta = theta(:);
phi = phi(:);