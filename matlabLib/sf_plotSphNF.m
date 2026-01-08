% Plots the the near field distribuition (psi) on a theta/phi plane
%
% sf_plotSphNF(matrixSize, theta, phi, psi, chosenTitle, dftFlag)
%
% IN: matrixSize = to collect the near field vector into a matrix with
%                  Ntheta rows and Nphi columns
%     theta-phi = vectors of the angles values
%     psi = vector of the near field sampled values
% INopt: chosenTitle = to put next to the |Psi| title
%        dftFlap = boolean for DFT coefficients plot if asserted
%
% Laurent Ntibarikure
function sf_plotSphNF(matrixSize, theta, phi, psi, chosenTitle, dftFlag)

if ~(nargin>4)
  chosenTitle = '';
  dftFlag = false;
elseif ~(nargin>5)
  dftFlag = false;
end
  

theta = vector2matrix(matrixSize, theta);
phi = vector2matrix(matrixSize, phi);
psi = vector2matrix(matrixSize, psi);

figure;
if ~dftFlag
  surf(theta*(180/pi), phi*180/pi, abs(psi),...
    'FaceAlpha',1,'EdgeAlpha',.1, 'EdgeColor','k');
  xlabel('\theta [°]');
  ylabel('\phi [°]');
  zlabel('|\Psi|')
  title(['|\Psi|', chosenTitle]);
else
  dftPsi = fftshift(fft2(psi));
  surf(theta*(180/pi), phi*180/pi, abs(dftPsi),...
    'FaceAlpha',.1,'EdgeAlpha',1, 'EdgeColor','k');
  xlabel('\theta [°]');
  ylabel('\phi [°]');
  zlabel('dft(|\Psi|)')
  title(['dft(|\Psi|)', chosenTitle]);
end
axis tight;
view([45 45]);