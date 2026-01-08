% Returns the box dimensions to encompass the array dimensions
%
% [xMin, xMax, yMin, yMax, zMin, zMax, xPts, yPts, zPts] = ...
%   getBoxDim(lambda, arrayPos, WLranging, WLspacing, ext )
%
% IN: lambda = wavelength
%     arrayPos = array elements positions of a planar array on XY plane
%     WLranging = distance between the array and the top and bottom faces
%                 of the enclosing box in wavelengths
%     WLspacing = sampling resolution on the box faces in wavelengths
%     ext = distance between the edge array elements and the lateral faces
%           in wavelengths
%
% OUT: xMin,xMax,yMin,yMax,zMin,zMax = dimensions of the box
%      xPts,yPts,zPts = number of sampling points along, respectively, 
%                       x, y and z directions
%
% Laurent Ntibarikure
function [xMin, xMax, yMin, yMax, zMin, zMax, xPts, yPts, zPts] = ...
  getBoxDim(lambda, arrayPos, WLranging, WLspacing, ext )

spacing = WLspacing * lambda;
ranging = WLranging * lambda;

xMin = min(arrayPos(1,:)) - ext * lambda;
xMax = max(arrayPos(1,:)) + ext * lambda;

yMin = min(arrayPos(2,:)) - ext * lambda;
yMax = max(arrayPos(2,:)) + ext * lambda;

zMin = min(arrayPos(3,:)) - ranging;
zMax = max(arrayPos(3,:)) + ranging;


xPts = floor( ( xMax-xMin ) / spacing);
yPts = floor( ( yMax-yMin ) / spacing);
zPts = floor( 2 * ranging / spacing);