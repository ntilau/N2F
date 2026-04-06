% Computes far-field directly from array sources (reference method)
%
% Direct summation for validation. Useful as ground truth for comparing
% with N2F transformation methods - both should give identical results.
%
% ALGORITHM: Far-field point source superposition
%   Theta-component: E_theta = sum_k excit_k * exp(ik*R_k) / (4*pi*R_k)
%   where R_k = source distance in far-field (plane wave approx)
%
% fPsiRef = sf_directffSolver(lambda, theta, phi, excitPhasor, arrayPos)
%
% IN: lambda = wavelength [m]
%     theta = polar far-field observation angles [rad] (nTheta x nPhi)
%     phi = azimuthal far-field observation angles [rad] (nTheta x nPhi)
%     excitPhasor = complex excitation amplitudes for each source (1 x nSources)
%     arrayPos = Cartesian source positions (3 x nSources)
%
% OUT: fPsiRef = far-field pattern (nPhi x nTheta)
%                Direct from sources - highest accuracy reference
%
% NOTE: O(nSources * nTheta * nPhi) complexity; use for validation/small problems
%       For large problems, use N2F operators (computed once, reused many times)
%
% Laurent Ntibarikure
function fPsiRef = sf_directffSolver(lambda, theta, phi, ...
  excitPhasor, arrayPos)
fprintf('#> Computing direct far fields ... ');
tic

kWN = 2*pi/lambda; % wavenumber
%% Reference patterns
fPsiRef = zeros(size(phi,2),size(theta,2));
for i=1:length(phi)
  fPsiRef(i,:) = excitPhasor / (4*pi) * exp(1i*kWN * ... 
    (arrayPos(1,:).' * sin(theta) * cos(phi(i)) + ...
    arrayPos(2,:).' * sin(theta) * sin(phi(i)) + ...
    arrayPos(3,:).' * cos(theta) ));
end

fprintf('%2.4g s.\n',toc);