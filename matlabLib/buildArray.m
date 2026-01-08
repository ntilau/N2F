% Returns the position of the radiating elements of a planar array built on
% the XY plane symmetrically centered in (x=0,y=0,z=0)
%
% arrayPos = buildArray(lambda, nbrElems_x, WLspacing_x, ...
%   nbrElems_y, WLspacing_y);
%
% IN: lambda = the wavelength in m
%     nbrElems_x = nbr of elements along the x direction
%     WLspacing_x = spacing in wavelengths between the elements along 
%                   the x direction
%     nbrElems_y = nbr of elements along the y direction
%     WLspacing_y = spacing in wavelengths between the elements along 
%                   the y direction
% OUT: arrayPos = array of the positions, every column correspond to x, y,
%                 z Cartesian coordinates [x;y;z]
%
% Example:
%   arrayPos = buildArray(.01, 4, .5, 2, .7) for the positions of a planar
%   array of 4 by 2 radiating elements equally spaced of 0.5*lambda along 
%   the x direction and 0.7*lambda along the y direction. lambda is of 1 cm
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