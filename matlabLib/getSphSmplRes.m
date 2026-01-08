% Computes the resolution in degrees given the desired resolution in
% wavelengths
%
% [dTheta, dPhi] = getSphSmplRes(radius, sphSmplRes)
%
% IN: radius = radius of the sphere in metres
%     sphSmplRes = resulotion in metres
%
% OUT: dTheta, dPhi = resolution in degres
%
% Laurent Ntibarikure
function [dTheta, dPhi] = getSphSmplRes(radius, sphSmplRes)

nTheta = floor(radius*pi/sphSmplRes-1);
nPhi = floor(radius*2*pi/sphSmplRes);

dTheta = 180/nTheta;
dPhi = 360/nPhi;

fprintf(' (dThetaNF=%2.4g°, dPhiNF=%2.4g°) ... ',dTheta, dPhi);