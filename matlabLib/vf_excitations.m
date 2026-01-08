% Computes the required linear phase excitations for a beam steering 
% process in (theta,phi). The phasors are computed in order to achieve
% constructive interference in the scan angle
%
% [J, M] =  vf_excitations(k0, arrayPos, Jmag, tJ, pJ,...
%   Mmag, tM, pM, steering_t, steering_p)
%
% IN: k0 = free-space wavenumber [1/m]
%     arrayPos = Cartesian coordinates of the array elements [m]
%     Jmag = electric current density
%     tJ-pJ = electric source orientation [°] (spherical components)
%     Mmag = magnetic current density
%     tM-pM = magnetic source orientation [°](spherical components)
%     steering_t = steering theta [°] (spherical coordinates)
%     steering_p = steering phi [°] (spherical coordinates)
%
% OUT: excitPhasor = vector of the excitations for unit power radiation and
%                    normalized with respect to a single isotropic radiator
%
% Laurent Ntibarikure
function [J, M] =  vf_excitations(k0, arrayPos, Jmag, tJ, pJ,...
  Mmag, tM, pM, steering_t, steering_p)

excitPhasor = exp(-1i*k0* ...
  (arrayPos(1,:).* sin(deg2rad(steering_t)) .* cos(deg2rad(steering_p)) + ...
  arrayPos(2,:).* sin(deg2rad(steering_t)) .* sin(deg2rad(steering_p)) + ...
  arrayPos(3,:).* cos(deg2rad(steering_t)) ));

[J(1,:),J(2,:),J(3,:)] = ...
  spherical2cartesian(excitPhasor*Jmag,0,0,deg2rad(tJ),deg2rad(pJ));
[M(1,:),M(2,:),M(3,:)] = ...
  spherical2cartesian(excitPhasor*Mmag,0,0,deg2rad(tM),deg2rad(pM));
