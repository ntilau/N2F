% Computes the required linear phase excitations for a beam steering 
% process in (theta,phi). The phasors are normalized such that the total
% power emitted is unitary
%
% excitPhasor = sf_excitations(lambda, arrayPos, steering_t, steering_p)
%
% IN: lambda = wavelength [m]
%     arrayPos = Cartesian coordinates of the array elements [m]
%     steering_t = steering theta [°] (spherical coordinates)
%     steering_p = steering phi [°] (spherical coordinates)
%
% OUT: excitPhasor = vector of the excitations for unit power radiation and
%                    normalized with respect to a single isotropic radiator
%
% Laurent Ntibarikure
function excitPhasor = sf_excitations(lambda, arrayPos, ...
  steering_t, steering_p)

fprintf('#> Computing excitations ... ');
tic

kWN = 2*pi/lambda; % wavenumber

excitPhasor = 4*pi/sqrt(size(arrayPos,2))* exp(-1i*kWN* ...
  (arrayPos(1,:).* sin(deg2rad(steering_t)) .* cos(deg2rad(steering_p)) + ...
  arrayPos(2,:).* sin(deg2rad(steering_t)) .* sin(deg2rad(steering_p)) + ...
  arrayPos(3,:).* cos(deg2rad(steering_t)) ));

fprintf('%2.4g s.\n',toc);