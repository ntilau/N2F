% Plots the near field solid (psi)
%
% sf_plotSphNFSolid(matrixSize, theta, phi, psi, chosenTitle)
%
% IN: matrixSize = to collect the near field vector into a matrix with
%                  Ntheta rows and Nphi columns
%     theta-phi = vectors of the angles values
%     psi = vector of the near field sampled values
% INopt: chosenTitle = to put next to the |Psi| title
%
% Laurent Ntibarikure
function sf_plotSphNFSolid(matrixSize, theta, phi, psi, chosenTitle)

offset=-2;
theta = vector2matrix(matrixSize, theta);
phi = vector2matrix(matrixSize, phi);
psi = vector2matrix(matrixSize, psi);

maxPsi=max(max(abs(psi)));
r=(1.1*abs(psi)/maxPsi);
r(r(:,:)<0)=0;

% [x, y, z] = sph2cart(phi,pi/2-theta, abs(psi)./ ...
%     max(max(abs(psi))));
[x, y, z] = sph2cart(phi,pi/2-theta, r);

figure;
clf();
surf(x, y, z, sqrt(x.^2+y.^2+z.^2),...
  'FaceAlpha',1,'EdgeAlpha',.3, 'EdgeColor','k');
xlabel('x');
ylabel('y');
zlabel('z')
axis([-1 1 -1 1 -1 1]);
axis('equal');
view([155 30]);
axis('off');
line([0 0],[0 0],[0 3+offset],'color','k');
line([0 0],[0 3+offset],[0 0],'color','k');
line([0 3+offset],[0 0],[0 0],'color','k');
text(0,0,3+offset,'z');
text(0,3+offset,0,'y');
text(3+offset,0,0,'x');            
title(['|\Psi| solid', chosenTitle]);