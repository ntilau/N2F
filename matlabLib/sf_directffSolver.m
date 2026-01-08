% Computes the far field by direct integration of the sources
%
% fPsiRef = sf_directffSolver(lambda, theta, phi, ...
%   excitPhasor, arrayPos)
%
% IN: lambda = wavelength
%     theta-phi = far field look angle
%     excitPhasor = phasors of the point sources
%     arrayPos = Cartesian coordinates of the point sources
%
% OUT: fPsiRef = far field detector for the look angles selected
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