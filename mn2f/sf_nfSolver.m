% Computes the scalar near field and its normal derivative over a
% surface used in a near-field to far-field transformation.
% The implementation uses point-source superposition and scalar Huygens'
% formulation.
%
% [psi, delPsi] = sf_nfSolver(lambda, excitPhasor, Rmag, NdotRV)
%
% IN:
%   lambda      = wavelength [m]
%   excitPhasor = 1xNs vector of complex source excitation phasors
%   Rmag        = Ns x Npts matrix of distances from each source to each
%                 surface sampling point
%   NdotRV      = Ns x Npts matrix of normal-direction dot products used
%                 for the normal derivative computation
%
% OUT:
%   psi   = 1xNpts near-field values on the surface
%   delPsi= 1xNpts normal derivative of the near field on the surface
%
% Note: psi and delPsi are computed using the free-space Green's function
%       for scalar point sources.
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