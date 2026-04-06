% Returns the positions of a planar radiating array defined on the XY
% plane, centered at the origin.
%
% arrayPos = buildArray(lambda, nbrElems_x, WLspacing_x, ...
%   nbrElems_y, WLspacing_y);
%
% IN:
%   lambda      = wavelength [m]
%   nbrElems_x  = number of elements along x
%   WLspacing_x = element spacing in wavelengths along x
%   nbrElems_y  = number of elements along y
%   WLspacing_y = element spacing in wavelengths along y
%
% OUT:
%   arrayPos = 3xN matrix of element coordinates [x; y; z], with z=0
%              and N = nbrElems_x * nbrElems_y.
%
% Example:
%   arrayPos = buildArray(.01, 4, .5, 2, .7)
%   creates a 4x2 planar array with 0.5 lambda spacing along x and
%   0.7 lambda spacing along y. The array is centered at the origin.
%
% Laurent Ntibarikure
function arrayPos = buildArray(lambda, nbrElems_x, WLspacing_x, ...
  nbrElems_y, WLspacing_y)

fprintf('#> Building array ... ');
tic

spacing_x = WLspacing_x*lambda;
spacing_y = WLspacing_y*lambda;
% checks odd or even for array centering
if mod(nbrElems_x,2)
  dim_x = floor(nbrElems_x/2);
else
  dim_x = (nbrElems_x-1)/2;
end
if mod(nbrElems_y,2)
  dim_y = floor(nbrElems_y/2);
else
  dim_y = (nbrElems_y-1)/2;
end
% create planar grid points for array elements positions
[x0,y0] = meshgrid(spacing_x*(-dim_x:dim_x), spacing_y*(-dim_y:dim_y));
arrayPos(1,:) = x0(:);
arrayPos(2,:) = y0(:);
arrayPos(3,:) = zeros(size(x0(:)));

fprintf('%2.4g s.\n', toc);