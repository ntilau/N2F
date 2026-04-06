% Computes far-field electric field components from near-field surface
% E and H fields using Huygens' principle for vector fields.
%
% [EtFF, EpFF] = vf_nf2ffSolver(k0, z0, surfPos, n, dS, E, H, thetaFF, phiFF)
%
% IN:
%   k0       = free-space wavenumber [rad/m]
%   z0       = free-space impedance [ohm]
%   surfPos  = 3xNpts Cartesian coordinates of the surface points [m]
%   n        = 3xNpts outward unit normals at each surface point
%   dS       = 1xNpts surface patch areas [m^2]
%   E        = 3xNpts electric near-field values [V/m]
%   H        = 3xNpts magnetic near-field values [A/m]
%   thetaFF  = Nphi x Ntheta (meshgrid) elevation angles [rad]
%   phiFF    = Nphi x Ntheta (meshgrid) azimuth angles [rad]
%
% OUT:
%   EtFF = Nphi x Ntheta theta-component of the far electric field [V/m]
%   EpFF = Nphi x Ntheta phi-component of the far electric field [V/m]
%
% This function computes equivalent surface electric and magnetic currents:
%   Js = n x H and Ms = E x n,
% then integrates their radiation contributions for each far-field direction.
%
% Laurent Ntibarikure
function [EtFF, EpFF] =  vf_nf2ffSolver(k0, z0, surfPos, n, dS,...
  E, H, thetaFF, phiFF)

tic();

% Equivalent surface currents for Huygens sources
Js = cross(n,H);
Ms = cross(E, n);

% Initialize far-field Cartesian components
ExFF = zeros(size(thetaFF));
EyFF = ExFF;
EzFF = ExFF;

for i = 1:size(thetaFF,2)
  for j = 1:size(phiFF,1)
    % Observation unit vector in the far-field direction
    RffVx = sin(thetaFF(j,i)).*cos(phiFF(j,i));
    RffVy = sin(thetaFF(j,i)).*sin(phiFF(j,i));
    RffVz = cos(thetaFF(j,i));

    % Phase factor for each surface point in the observation direction
    green = exp(1i*k0*(RffVx.*surfPos(1,:) + ...
      RffVy.*surfPos(2,:) + RffVz.*surfPos(3,:)));

    % Projection of equivalent electric current Js on the observation direction
    JsdotRffV = Js(1,:).*RffVx + Js(2,:).*RffVy + Js(3,:).*RffVz;

    [MscrossRffVx,MscrossRffVy,MscrossRffVz]=...
      crossOperator(Ms(1,:),Ms(2,:),Ms(3,:),RffVx,RffVy,RffVz);

    % Surface integral for the far-field Cartesian components
    ExFF(j,i)= ExFF(j,i) - 1i*k0/(4*pi).*sum( ...
      (z0*(Js(1,:)-JsdotRffV.*RffVx) + MscrossRffVx).*green.*dS );
    EyFF(j,i)= EyFF(j,i) - 1i*k0/(4*pi).*sum( ...
      (z0*(Js(2,:)-JsdotRffV.*RffVy) + MscrossRffVy).*green.*dS );
    EzFF(j,i)= EzFF(j,i) - 1i*k0/(4*pi).*sum( ...
      (z0*(Js(3,:)-JsdotRffV.*RffVz) + MscrossRffVz).*green.*dS );
  end
end

fprintf('#> nf2ff computation time : %g s.\n',toc)

% Convert Cartesian far fields to spherical components. The algorithm
% assumes the standard far-field approximation and omits the radial 1/R
% decay factor for relative pattern evaluation.
[ErFF,EtFF,EpFF] = cartesian2spherical(ExFF, EyFF, EzFF, thetaFF, phiFF);