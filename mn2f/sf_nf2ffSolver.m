% Computes the scalar near-field to far-field transformation using
% Huygens' principle over a sampled bounding surface.
%
% fPsi = sf_nf2ffSolver(lambda, theta, phi, surfPos, N, dS, psi, delPsi)
%
% IN:
%   lambda  = wavelength [m]
%   theta    = 1xNtheta vector of far-field elevation angles [rad]
%   phi      = 1xNphi vector of far-field azimuth angles [rad]
%   surfPos  = 3xNpts coordinates of surface sampling points [m]
%   N        = 3xNpts outward unit normal vectors at each surface point
%   dS       = 1xNpts surface patch areas
%   psi      = 1xNpts scalar near-field values at surface points
%   delPsi   = 1xNpts scalar normal derivative values at surface points
%
% OUT:
%   fPsi = Nphi x Ntheta far-field pattern matrix for specified angles
%
% The algorithm computes the surface integral of the Huygens sources
% for each observation direction (theta,phi).
%
% Laurent Ntibarikure
function fPsi = sf_nf2ffSolver(lambda, theta, phi, ...
  surfPos, N, dS, psi, delPsi)
fprintf('#> Computing n2f fields transformations ... ');
tic

kWN = 2*pi/lambda; % wavenumber
%% Computing far fields by numerical integration on the sampled surface
fPsi = zeros(size(phi,2), size(theta,2));
for i=1:length(theta)
  for j=1:length(phi)
    % Unit direction vector for the far-field observation
    RxV = sin(theta(i)) .* cos(phi(j));
    RyV = sin(theta(i)) .* sin(phi(j));
    RzV = cos(theta(i));
    R = [RxV, RyV, RzV];

    % Projection of surface normals onto the far-field direction
    nR(1,:) = RxV .* N(1,:);
    nR(2,:) = RyV .* N(2,:);
    nR(3,:) = RzV .* N(3,:);
    ndotR = sum(nR,1);

    % Scalar Green's function for each surface point in the far-field direction
    green = 1/(4*pi)*exp(1i*kWN.*(R*surfPos));

    % Surface integral for the far-field pattern component
    fPsi(j,i) = sum( green.*( (1i*kWN.*ndotR).*psi - delPsi) .* dS);
  end
end

fprintf('%2.4g s.\n',toc);