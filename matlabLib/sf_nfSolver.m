% Computes the scalar near field and its derivative over the surface
% chosen for the near field to far field computation (scalar Huygens' 
% principle)
%
% [psi, delPsi]=  sf_nfSolver(lambda, excitPhasor, Rmag, NdotRV)
%
% IN: lambda = wavelength [m]
%     excitPhasor = phasors of the excitations of the point sources
%     Rmag = for each pt source, distances to the sampling points on the
%            bounding surface
%     NdotRV = dot product for derivatives computation
%
% OUT: psi = near field (normalized for pattern computation by the chosen 
%            sources)
%      delPsi = near field normal derivative (also normalized)
%
% Laurent Ntibarikure
function [psi, delPsi]=  sf_nfSolver(lambda, excitPhasor, Rmag, NdotRV)

fprintf('#> Computing near fields and derivatives ... ');
tic
kWN = 2*pi/lambda; % wavenumber
%% Computing psi, the near field
psi = excitPhasor * (exp(-1i*kWN*Rmag)./(4*pi*Rmag));
%% Computing delPsi, the near field normal derivative
delPsi = excitPhasor * (-(1i*kWN+1./Rmag).* ...
  exp(-1i*kWN*Rmag)./(4*pi*Rmag).*NdotRV);
fprintf('%2.4g s.\n',toc);