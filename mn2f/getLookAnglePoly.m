% Returns the complex trigonometric polynomials values for selected look
% angles and related to the proper eigenfrequency
%
% [ifftOp, phiFF, thetaFF] = getLookAnglePoly(phi, sizet, ...
%   sizetref, nbrCoeffs)
%
% IN: phi = angles of cut planes selected
%     sizet = number of theta look angles for pattern plot
%     sizetref = size of the initial DFT dimension
%     nbrCoeffs = number of coefficients retained in DFT-truncation
%
% OUT: ifftOp = complex trigonometric polynomials sampled values matrix
%      phiFF, thetaFF = sampled look angles with phi and sizet params
%
% Laurent Ntibarikure
function [ifftOp, phiFF, thetaFF] = getLookAnglePoly(phi, sizet, ...
  sizetref, nbrCoeffs)

nbrCoeffs = floor(nbrCoeffs/2);
[phiFF,thetaFF]=meshgrid(phi, ...
   linspace(0,2*pi*(sizet)/sizet,sizet));
angles = linspace(0,(sizet-1)/sizet,sizet);
ifftOp=zeros(sizet,2*nbrCoeffs+1);
for i=1:nbrCoeffs+1
    ifftOp(:,i)= sizetref\exp(1i*2*pi*angles*(i-1));
end
for i=1:nbrCoeffs
   ifftOp(:,nbrCoeffs+1+i)= sizetref\exp(-1i*2*pi*angles*(nbrCoeffs-i+1));
end