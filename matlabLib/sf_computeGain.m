% Returns the directivity in [dBi] for properly normalized far field
%
% gain = sf_computeGain(fPsi)
%
% IN: psi = power density normalized far field detectors
%
% OUT: gain = directivity in [dBi]
%
% Laurent Ntibarikure
function gain = sf_computeGain(fPsi)

gain = abs(fPsi).^2;