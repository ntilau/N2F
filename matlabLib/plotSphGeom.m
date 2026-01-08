% Plots the sphere and array geometry
%
% plotSphGeom(matrixSize, arrayPos, spherePos)
%
% IN: matrixSize = for theta-phi assembly
%     arrayPos = for plot of the array positions
%     spherePos = for plot of the sphere grid (surface patches) n.b.
%                 sampling occurs in the middle of the patches
%
% Laurent Ntibarikure
function plotSphGeom(matrixSize, arrayPos, spherePos)

X = vector2matrix(matrixSize,spherePos(1,:).');
Y = vector2matrix(matrixSize,spherePos(2,:).');
Z = vector2matrix(matrixSize,spherePos(3,:).');

figure;
colormap(getColorMap);
mesh(X, Y, Z, 'FaceAlpha',0,'EdgeAlpha',1, 'EdgeColor', getColorMap);
% grid off;
hidden off
hidden on
hold on;
xlabel('x [\lambda]');
ylabel('y [\lambda]');
zlabel('z [\lambda]')
title('Sphere geometry');
axis tight;
axis equal;
view(45,45);
for i=1:length(arrayPos)
  plot3(arrayPos(1,i), arrayPos(2,i), arrayPos(3,i), '.k', 'MarkerSize',20);
  hold on;
end