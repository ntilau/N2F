% Constructs a spherical bounding surface with flexible angular sampling strategies
%
% Generates regularly or adaptively sampled discrete sphere for N2F transformations.
% Key for efficient far-field computation via spherical Huygens' principle.
%
% SAMPLING MODES (flag parameter):
%   flag = 0 : Standard dTheta-dPhi uniform grid (most common)
%              Regular Cartesian grid mapped to sphere
%   flag = 1 : Enhanced resolution using wavelength-based sampling (sphSmplRes)
%              Automatically adjusts density based on electrical size
%   flag = 2 : Special dTheta-dPhi grid optimized for 3D visualization
%              Avoids patches overlapping in plot (esthetic improvement)
%              WARNING: May introduce small N2F errors; use flag=0 for accurate computation
%   flag = 3 : Custom adaptive sampling (user-defined via rotAngle, rotAxis)
%              Applies rotation for directional focus or regional refinement
%
% [spherePos, dS, theta, phi, matrixSize] = ...
%   buildSphere(radius, sphSmplRes, dTheta, dPhi, flag, rotAngle, rotAxis)
%
% IN: radius = sphere radius [m]
%     sphSmplRes = wavelength-based sampling resolution [wavelengths]
%                  (used if flag=1; ignored otherwise)
%     dTheta = polar angle discretization [degrees] (typical: 1-5 degrees)
%     dPhi = azimuthal angle discretization [degrees] (typical: 1-5 degrees)
%     flag = sampling mode selector (0,1,2,3)
%     rotAngle = rotation angle [degrees] for flag=3 mode (e.g., 45 for tilt)
%     rotAxis = unit vector [x;y;z] defining rotation axis (e.g., [1;1;0] for diagonal)
%
% OUT: spherePos = Cartesian coordinates on sphere surface (3 x nNodes)
%                  All points satisfy ||spherePos|| = radius identically
%      dS = area per sampling patch [m^2] (1 x nNodes)
%           For uniform sphere: approx (4*pi*radius^2) / nNodes
%      theta = polar angles [radians] of sampling nodes (0 to pi)
%      phi = azimuthal angles [radians] of sampling nodes (0 to 2*pi)
%      matrixSize = mesh metadata: [nTheta x nPhi], total points, etc.
%
% EFFICIENCY: Sphere avoids general poly geometry distance computations.
%   All properties (point locations, normals, areas) have closed-form expressions,
%   enabling O(1) per-point computation vs. O(n) for arbitrary surfaces.
%
% Laurent Ntibarikure
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
