% Computes normalized directivity gain from complex far-field amplitude
%
% DEFINITION: Gain = 4*pi*|E|^2 / P_total  [dimensionless]
%             In dB: Gain_dB = 10*log10(|E|^2) - normalized peak
%
% Used for pattern comparison, directivity assessment, and error metrics.
% Converts single polarization component to power-based gain figure.
%
% gain_dB = sf_computeGain(fPsi)
%
% IN: fPsi = complex far-field amplitude (output from N2F or direct solver)
%            Single polarization component (theta or phi)
%
% OUT: gain_dB = normalized directivity gain [dB]
%                Max value = 0 dB (at main lobe peak)
%                Sidelobe levels relative to main lobe
%
% Laurent Ntibarikure
function gain = sf_computeGain(fPsi)

gain = abs(fPsi).^2;