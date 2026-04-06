function[Ax,Ay,Az] = spherical2cartesian(Ar,At,Ap,theta,phi)
Ax = Ar.*sin(theta).*cos(phi) + At.*cos(theta).*cos(phi) - ...
    Ap.*sin(phi);
Ay = Ar.*sin(theta).*sin(phi) + At.*cos(theta).*sin(phi) + ...
    Ap.*cos(phi);
    Az = Ar.*cos(theta) - At.*sin(theta);
endfunction