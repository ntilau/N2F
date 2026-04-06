% Returns the gain(P=Pacc)-directivity(P=Prad)in theta and phi polarizations
%
% [gaint, gainp] = vf_computeGain(z0, EtFF, EpFF, P)
%
% IN: z0 = free-space wave impedance
%     EtFF, EpFF = far electric field detector in spherical coordinates
%     P = power for comparison
%
% OUT: gaint, gainp = gain related to theta and phi polarizations
%                    total gain = gaint + gainp
%
% Laurent Ntibarikure
function [gaint, gainp] = vf_computeGain(z0, EtFF, EpFF, P)

gaint=2*pi/P.*abs(EtFF).^2./z0; % gain theta polarization
gainp=2*pi/P.*abs(EpFF).^2./z0; % gain phi polarization