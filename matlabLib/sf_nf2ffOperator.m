% Computes scalar near-field-to-far-field operator matrices using
% Huygens' principle.
%
% [Lpsi, LdelPsi] = sf_nf2ffOperator(lambda, theta, phi, surfPos, N, dS)
%
% IN:
%   lambda  = wavelength [m]
%   theta   = 1xNtheta vector of far-field elevation angles [rad]
%   phi     = 1xNphi vector of far-field azimuth angles [rad]
%   surfPos = 3xNpts Cartesian coordinates of bounding surface points
%   N       = 3xNpts outward normal unit vectors at each point
%   dS      = 1xNpts surface patch areas
%
% OUT:
%   Lpsi    = Ntheta x Npts x Nphi operator for scalar field values
%   LdelPsi = Ntheta x Npts x Nphi operator for normal derivatives
%
% These operators allow the far-field pattern to be computed as:
%   fPsi(:,i) = Lpsi(:,:,i) * psi.' + LdelPsi(:,:,i) * delPsi.'
% for each azimuth slice phi(i).
%
% Laurent Ntibarikure
function [Lpsi, LdelPsi] = sf_nf2ffOperator(lambda, theta, phi, ...
  surfPos, N, dS)
fprintf('#> Computing n2f operators ... ');
tic
kWN = 2*pi/lambda; % wavenumber
%% Build operator matrices for each phi slice
Lpsi = zeros(size(theta,2), size(surfPos,2), size(phi,2));
LdelPsi = Lpsi;
for i=1:length(phi)
  % Direction cosines for the far-field observation directions
  RxV = sin(theta.') .* cos(phi(i));
  RyV = sin(theta.') .* sin(phi(i));
  RzV = cos(theta.');
  R = [RxV, RyV, RzV];

  % Green's function times surface element area for each direction and point
  greendS = 1/(4*pi)*exp(1i*kWN.*R*surfPos) .* (ones(size(theta,2),1) * dS);

  % Operator components for psi and normal derivative contributions
  Lpsi(:,:,i) = greendS .* (1i*kWN.*(R*N));
  LdelPsi(:,:,i) = - greendS;
end

fprintf('%2.4g s.\n',toc);
