% trasformation from spherical components to cartesian components of
% arbitrary structured vectorial field A
%
% [Ax,Ay,Az] = spherical2cartesian(Ar,At,Ap,theta,phi)
%
% IN: Ar, At, Ap = spherical components of A
%     theta, phi = look angle (scalar)
%
% OUT: Ax, Ay, Az = Cartesian components of A
%
% Laurent Ntibarikure
function [Ax,Ay,Az] = spherical2cartesian(Ar,At,Ap,theta,phi)
Ax = Ar.*sin(theta).*cos(phi) + At.*cos(theta).*cos(phi) - ...
    Ap.*sin(phi);
Ay = Ar.*sin(theta).*sin(phi) + At.*cos(theta).*sin(phi) + ...
    Ap.*cos(phi);
Az = Ar.*cos(theta) - At.*sin(theta);