% Computes the near field to far field transformation using scalar Huygens'
% principle
%
% fPsi = sf_nf2ffSolver(lambda, theta, phi, ...
%   surfPos, N, dS, psi, delPsi)
%
% IN: lambda = wavelength
%     theta-phi = far field look angle
%     surfPos = bounding surface sampling points Cartesian coordinates
%     N = outwardly directed normal unit vector for each surface sampling
%         point
%     dS = surface patches area
%     psi-delPsi = near field and normal derivative on the sampling points
%
% OUT: fPsi = far field detector for the look angles selected
%
% Laurent Ntibarikure
function fPsi = sf_nf2ffSolver(lambda, theta, phi, ...
  surfPos, N, dS, psi, delPsi)
fprintf('#> Computing n2f fields transformations ... ');
tic

kWN = 2*pi/lambda; % wavenumber
%% Computing far fields
fPsi = zeros(size(phi,2), size(theta,2));
for i=1:length(theta)
  for j=1:length(phi)
    RxV = sin(theta(i)) .* cos(phi(j));
    RyV = sin(theta(i)) .* sin(phi(j));
    RzV = cos(theta(i));
    R = [RxV, RyV, RzV];
    nR(1,:) = RxV .* N(1,:);
    nR(2,:) = RyV .* N(2,:);
    nR(3,:) = RzV .* N(3,:);
    ndotR = sum(nR,1);
    green = 1/(4*pi)*exp(1i*kWN.*(R*surfPos));
    fPsi(j,i) = sum( green.*( (1i*kWN.*(ndotR)).*psi - delPsi) .* dS);
  end
end

fprintf('%2.4g s.\n',toc);