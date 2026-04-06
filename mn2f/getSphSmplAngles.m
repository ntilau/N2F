% Computes the sphere sampling angles for near field to far field
% transformation
%
% [theta, phi, matrixSize] = getSphSmplAngles(dTheta, dPhi)
%
% IN: dTheta, dPhi = theta and phi resolutions in degrees
% 
% OUT: theta, phi = column vectors of the angles
%      matrixSize = size of the angles meshgrid for plots
%
% Laurent Ntibarikure
function [theta, phi, matrixSize] = getSphSmplAngles(dTheta, dPhi)

[theta, phi] = meshgrid(deg2rad(dTheta/2:dTheta:180-dTheta/2),...
  deg2rad(0:dPhi:360-dPhi/2));
matrixSize = size(theta);
theta = theta(:);
phi = phi(:);