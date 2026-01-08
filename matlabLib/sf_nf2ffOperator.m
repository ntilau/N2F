% Computes the near field to far field transformation operator based on 
% scalar Huygens' principle
%
% [Lpsi, LdelPsi] = sf_nf2ffOperator(lambda, theta, phi, ...
%   surfPos, N, dS)
%
% IN: lambda = wavelength
%     theta-phi = far field look angle
%     surfPos = bounding surface sampling points Cartesian coordinates
%     N = outwardly directed normal unit vector for each surface sampling
%         point
%     dS = surface patches area
%
% OUT: Lpsi, LdelPsi = near field to far field operators for the look 
%                      angles selected
%
% Laurent Ntibarikure
function [Lpsi, LdelPsi] = sf_nf2ffOperator(lambda, theta, phi, ...
  surfPos, N, dS)
fprintf('#> Computing n2f operators ... ');
tic
kWN = 2*pi/lambda; % wavenumber
%% Computing  radiation pattern, the directivity (gain)
Lpsi = zeros(size(theta,2), size(surfPos,2), size(phi,2));
LdelPsi = Lpsi;
for i=1:length(phi)
  RxV = sin(theta.') .* cos(phi(i));
  RyV = sin(theta.') .* sin(phi(i));
  RzV = cos(theta.');
  R = [RxV, RyV, RzV];
  greendS = 1/(4*pi)*exp(1i*kWN.*R*surfPos) .* (ones(size(theta,2),1) * dS);
  Lpsi(:,:,i) = greendS .* (1i*kWN.*(R*N));
  LdelPsi(:,:,i) =  - greendS;
end

fprintf('%2.4g s.\n',toc);
