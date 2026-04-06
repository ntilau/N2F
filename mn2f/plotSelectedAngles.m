% plots the unit vectors pointing in the scan directions of both spanning
% vectors and testing vectors
%
% plotSelectedAngles(spanT, spanP, testT, testP)
%
% IN: spanT, spanP = theta and phi of the scanning direction
%     testT, testP = theta and phi of the testing direction
%
% Laurent Ntibarikure
function plotSelectedAngles(spanT, spanP, ...
  testT, testP)

figProp = getFigureProperties();
map = getColorMap();
x = sin(deg2rad(spanT)) .* cos(deg2rad(spanP));
y = sin(deg2rad(spanT)) .* sin(deg2rad(spanP));
z = cos(deg2rad(spanT));
o = zeros(size(x));
if nargin > 2
  xt = sin(deg2rad(testT)) .* cos(deg2rad(testP));
  yt = sin(deg2rad(testT)) .* sin(deg2rad(testP));
  zt = cos(deg2rad(testT));
  ot = zeros(size(xt));
end

figure; quiver3(o,o,o,x,y,z,1, 'LineWidth', figProp.lw, ...
  'MarkerSize',figProp.ms, 'Color', map, 'HandleVisibility', 'off');
hold on;
plot3(x,y,z,'.k', 'MarkerSize',25); axis('equal');
if nargin > 2
  hold on; quiver3(ot,ot,ot,xt,yt,zt,1, 'LineWidth', figProp.lw, ...
  'MarkerSize',figProp.ms, 'Color', 'r', 'HandleVisibility', 'off');
  hold on;
  plot3(xt,yt,zt,'sr', 'MarkerSize',10); axis('equal');
  legend('Scan angle space selection', 'Tested scan angles');
else
  legend('Scan angle space selection');
end
view(120,30);
axis('equal'); axis('tight');% axis('off')
grid off;
xlabel('x', 'FontSize', figProp.fs);
ylabel('y', 'FontSize', figProp.fs);
zlabel('z', 'FontSize', figProp.fs);