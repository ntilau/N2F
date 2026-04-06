% Computes the direct far-field pattern of a set of point sources.
% This reference solution is used to validate the near-field-based
% N2F transformation.
%
% fPsiRef = sf_directffSolver(lambda, theta, phi, excitPhasor, arrayPos)
%
% IN:
%   lambda       = wavelength [m]
%   theta, phi   = 1xNtheta and 1xNphi far-field angles [rad]
%   excitPhasor  = 1xNa complex excitation phasors for each source
%   arrayPos     = 3xNa Cartesian coordinates of the point sources [m]
%
% OUT:
%   fPsiRef = Nphi x Ntheta far-field pattern from the point source array
%
% The far-field pattern is computed by direct summation of the source
% contributions using the plane-wave phase factor for each observation
% direction.
%
% Laurent Ntibarikure
function fPsiRef = sf_directffSolver(lambda, theta, phi, ...
  excitPhasor, arrayPos)
fprintf('#> Computing direct far fields ... ');
tic

kWN = 2*pi/lambda; % wavenumber
%% Reference patterns
fPsiRef = zeros(size(phi,2),size(theta,2));
for i=1:length(phi)
  fPsiRef(i,:) = excitPhasor / (4*pi) * exp(1i*kWN * ... 
    (arrayPos(1,:).' * sin(theta) * cos(phi(i)) + ...
    arrayPos(2,:).' * sin(theta) * sin(phi(i)) + ...
    arrayPos(3,:).' * cos(theta) ));
end

fprintf('%2.4g s.\n',toc);