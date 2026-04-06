% Computes the near electric and magnetic fields produced by arrays of
% electric and magnetic currents over a bounding surface.
% The output is the surface E and H fields, Poynting flux, and total
% radiated power.
%
% [E, H, S, Pr] = vf_nfSolver(k0, z0, arrayPos, surfPos, n, dS, J, M)
%
% IN:
%   k0       = free-space wavenumber [rad/m]
%   z0       = free-space impedance [ohm]
%   arrayPos = 3xNa positions of the radiating current sources [m]
%   surfPos  = 3xNpts positions of the surface sampling points [m]
%   n        = 3xNpts outward unit normals at the surface points
%   dS       = 1xNpts surface patch areas [m^2]
%   J        = 3xNa electric current densities at each source
%   M        = 3xNa magnetic current densities at each source
%
% OUT:
%   E  = 3xNpts electric field on the surface [V/m]
%   H  = 3xNpts magnetic field on the surface [A/m]
%   S  = 1xNpts Poynting vector normal component on the surface
%   Pr = scalar total radiated power [W]
%
% The implementation is based on Kottler's equations for near-field
% source contributions from electric and magnetic currents.
%
% Laurent Ntibarikure
function [E, H, S, Pr] = vf_nfSolver(k0, z0, ...
  arrayPos, surfPos, n, dS, J, M)

% Initialize surface field arrays
E = zeros(size(surfPos));
H = zeros(size(surfPos));

for i=1:size(arrayPos,2)
  % Vector from the i-th source to each surface sampling point
  Rx = surfPos(1,:) - arrayPos(1,i);
  Ry = surfPos(2,:) - arrayPos(2,i);
  Rz = surfPos(3,:) - arrayPos(3,i);
  R = sqrt(Rx.^2 + Ry.^2 + Rz.^2);
  RxV = Rx ./ R; RyV = Ry ./ R; RzV = Rz ./ R;

  % Source current projections along the source-to-surface direction
  JdotRV = J(1,i).*RxV + J(2,i).*RyV + J(3,i).*RzV;
  MdotRV = M(1,i).*RxV + M(2,i).*RyV + M(3,i).*RzV;

  % Cross-products used in the Kottler field expressions
  [JcrossRVx,JcrossRVy,JcrossRVz]= ...
    crossOperator(J(1,i),J(2,i),J(3,i),RxV,RyV,RzV);
  [McrossRVx,McrossRVy,McrossRVz]= ...
    crossOperator(M(1,i),M(2,i),M(3,i),RxV,RyV,RzV);

  % Electric field contributions from electric and magnetic currents
  E(1,:) = E(1,:) + ...
    ( z0*k0^2/(4*pi).*( J(1,i).*(-1i./(k0*R)-1./(k0*R).^2+1i./(k0*R).^3) +...
    JdotRV.*RxV.*(1i./(k0*R)+3./(k0*R).^2-3i./(k0*R).^3) ) - k0^2/(4*pi).*...
    McrossRVx.*(1i./(k0.*R)+1./(k0.*R).^2) ).*exp(-1i*k0.*R);
  E(2,:) = E(2,:) + ...
    ( z0*k0^2/(4*pi).*( J(2,i).*(-1i./(k0*R)-1./(k0*R).^2+1i./(k0*R).^3) +...
    JdotRV.*RyV.*(1i./(k0*R)+3./(k0*R).^2-3i./(k0*R).^3) ) - k0^2/(4*pi).*...
    McrossRVy.*(1i./(k0.*R)+1./(k0.*R).^2) ).*exp(-1i*k0.*R);
  E(3,:) = E(3,:) + ...
    ( z0*k0^2/(4*pi).*( J(3,i).*(-1i./(k0*R)-1./(k0*R).^2+1i./(k0*R).^3) +...
    JdotRV.*RzV.*(1i./(k0*R)+3./(k0*R).^2-3i./(k0*R).^3) ) - k0^2/(4*pi).*...
    McrossRVz.*(1i./(k0.*R)+1./(k0.*R).^2) ).*exp(-1i*k0.*R);

  % Magnetic field contributions from magnetic and electric currents
  H(1,:) = H(1,:) + ...
    ( k0^2/(4*pi*z0).*( M(1,i).*(-1i./(k0*R)-1./(k0*R).^2+1i./(k0*R).^3) +...
    MdotRV.*RxV.*(1i./(k0*R)+3./(k0*R).^2-3i./(k0*R).^3) ) + k0^2/(4*pi).*...
    JcrossRVx.*(1i./(k0.*R)+1./(k0.*R).^2) ) .*exp(-1i*k0.*R);
  H(2,:) = H(2,:) + ...
    ( k0^2/(4*pi/z0).*( M(2,i).*(-1i./(k0*R)-1./(k0*R).^2+1i./(k0*R).^3) +...
    MdotRV.*RyV.*(1i./(k0*R)+3./(k0*R).^2-3i./(k0*R).^3) ) + k0^2/(4*pi).*...
    JcrossRVy.*(1i./(k0.*R)+1./(k0.*R).^2) ).*exp(-1i*k0.*R);
  H(3,:) = H(3,:) + ...
    ( k0^2/(4*pi*z0).*( M(3,i).*(-1i./(k0*R)-1./(k0*R).^2+1i./(k0*R).^3) +...
    MdotRV.*RzV.*(1i./(k0*R)+3./(k0*R).^2-3i./(k0*R).^3) ) + k0^2/(4*pi).*...
    McrossRVz.*(1i./(k0.*R)+1./(k0.*R).^2) ).*exp(-1i*k0.*R);
end

% Compute the normal component of the complex Poynting vector and total power
Sp = cross(E,conj(H));
S = dot(Sp,n);
Pr = 1/2 * real(sum(S .* dS));

fprintf('##> Radiated power : %g Watts.\n',Pr);