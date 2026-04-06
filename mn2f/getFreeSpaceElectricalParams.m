% Returns the electrical parameters given a specific frequency
%
% [z0,k0,lambda0] = getFreeSpaceElectricalParams(freq, unitary)
%
% IN: freq = frequency of analysis
%     unitary = returns values such that the wavelength is unitary
%
% OUT: z0 = free-space wave impedance = sqrt(mu0/eps0)
%      k0 = free-space wavenumber = 2*pi*freq*sqrt(eps0*mu0)
%      lambda0 = free-space wavelength = 2*pi/k0
%
% Laurent Ntibarikure
function [z0,k0,lambda0] = getFreeSpaceElectricalParams(freq, unitary )

eps0=8.85418781761e-12;
mu0=pi*4e-7;
c0=1/sqrt(eps0*mu0);
if nargin>1 && unitary
  freq = c0;
end
z0=sqrt(mu0/eps0);
lambda0=c0/freq;
k0=2*pi*freq/c0;