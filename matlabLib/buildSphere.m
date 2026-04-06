% Builds the near field sampling points location and the area of the
% patches
%
% [spherePos, dS, theta, phi, matrixSize] = ...
%   buildSphere(radius, sphSmplRes, dTheta, dPhi, flag, ...
%   rotAngle, rotAxis);
%
% IN: radius = radius of the sphere in [m]
%     sphSmplRes = sampling resolution in wavelengths (depend on flag)
%     dTheta = theta sampling resolution in degrees (depend on flag)
%     dPhi = phi sampling resolution in degrees (depend on flag)
%     flag = 0: dTheta-dPhi  1: sphSmplRes 2: dTheta-dPhi for nf plots
%            (chosen in a range that allows 3D plots without missing 
%            patches) n.b. 2: induces errors in the n2f computation
%     rotAngle = rotation angle in degrees of the sphere points for nf
%                rotation
%     rotAxis = rotation axis defined by the vector with components [x;y;z]
%
% OUT: spherePos = position [x;y;z] of the sphere sampling points
%      dS = area of the patch in which the fields are sampled
%      theta, phi = direction of the sphere sampling point
%      matrixSize = for the collection of the sampling points in terms of
%                   the sampling direction [theta, phi]
%
% Laurent Ntibarikure
function [spherePos, dS, theta, phi, matrixSize] = ...
  buildSphere(radius, sphSmplRes, dTheta, dPhi, flag, rotAngle, rotAxis)
fprintf('#> Building sphere ...');
tic

if flag == 1 || flag == 3
  [dTheta, dPhi] = getSphSmplRes(radius, sphSmplRes);
end
if flag == 2 || flag == 3% for plots (leading to erroneous surface integration)
  [theta, phi, matrixSize] = getSphSmplAnglesForPlots(dTheta, dPhi);
else
  [theta, phi, matrixSize] = getSphSmplAngles(dTheta, dPhi);
end
%% Sphere patches areas
dS(1,:) = radius.^2 .* sin(theta) * 2*pi^2 / size(theta,2) / size(phi,1);
%% Sphere patches vectors in cartesian coordinates
spherePos(1,:) = radius * sin(theta) .* cos(phi);
spherePos(2,:) = radius * sin(theta) .* sin(phi);
spherePos(3,:) = radius * cos(theta);
%% rotation
if nargin > 5
  rotMatrix = getRotationMatrix(rotAngle,rotAxis);
  spherePos = rotMatrix * spherePos;
end

fprintf('%2.4g s.\n',toc);
