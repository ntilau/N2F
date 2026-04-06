% trasformation from cartesian components to spherical components of
% arbitrary structured (scalar, vector, matrix) vectorial field A
%
% [Ar,At,Ap] = cartesian2spherical(Ax,Ay,Az,theta,phi)
%
% IN: Ax, Ay, Az = Cartesian components of A
%     theta, phi = look angle (scalar)
%
% OUT: Ar, At, Ap = spherical components of A
%
% Laurent Ntibarikure
function [Ar,At,Ap] = cartesian2spherical(Ax,Ay,Az,theta,phi)
Ar = Ax.*sin(theta).*cos(phi) + Ay.*sin(theta).*sin(phi) + ...
    Az.*cos(theta);
At = Ax.*cos(theta).*cos(phi) + Ay.*cos(theta).*sin(phi) - ...
    Az.*sin(theta);
Ap = - Ax.*sin(phi) + Ay.*cos(phi);