% Plots the array point sources
%
% sf_plotArrayGeom(arrayPos)
%
% IN: arrayPos = point sources positions in Cartesian coordinates
%
% Laurent Ntibarikure
function sf_plotArrayGeom(arrayPos)

figure(gcf);
hidden off
hidden on
hold on;
xlabel('x [\lambda]');
ylabel('y [\lambda]');
zlabel('z [\lambda]')
axis tight;
axis equal;
view(45,45);
for i=1:length(arrayPos)
  plot(arrayPos(1,i), arrayPos(2,i), '.k', 'MarkerSize',20);
  hold on;
end