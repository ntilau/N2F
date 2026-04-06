//% trasformation from cartesian components to spherical components of
//% arbitrary structured vectorial field A
function[Ar,At,Ap] = cartesian2spherical(Ax,Ay,Az,theta,phi)
Ar = Ax.*sin(theta).*cos(phi) + Ay.*sin(theta).*sin(phi) + ...
    Az.*cos(theta);
At = Ax.*cos(theta).*cos(phi) + Ay.*cos(theta).*sin(phi) - ...
    Az.*sin(theta);
Ap = - Ax.*sin(phi) + Ay.*cos(phi);
endfunction