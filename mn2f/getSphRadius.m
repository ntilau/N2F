% Computes the radius of the sphere provided the array elements positions 
% and the extention in wavelengths from the edge array element
%
% radius = getSphRadius(lambda, arrayPos, ext)
%
% IN: lambda = wavelength in [m]
%     arrayPos = coordinates of the array elements
%     ext = distance in wavelengths between the sphere and the edge array
%           elements
%
% OUT: radius = sphere radius in [m]
%
% Laurent Ntibarikure
function radius = getSphRadius(lambda, arrayPos, ext)

radius =  (sqrt(max(abs(arrayPos(1,:)))^2 + max(abs(arrayPos(2,:)))^2 + ...
  max(abs(arrayPos(3,:)))^2) / ...
  lambda + ext) * lambda;