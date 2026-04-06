% Spherical spiraling trajectory for septentrional hemisphere scanning
%
% [trajTheta, trajPhi, oppTrajTheta, oppTrajPhi] = ...
%   getSpiralingHelicoidalTrajectory(N, M, plot, fact)
%
% IN: N = number of scanning angles
%     M = number of test angles in the opposite trajectory
%     plot = if asserted plots the trajectory
%     fact = factor for theta range. if fact==1, theta is in the range of
%            [0°,90°], fact==1/2 -> [0°,45°]
%
% OUT: trajTheta, trajPhi = scan angles coordinates
%      oppTrajTheta, oppTrajPhi = test scan angles coordinates on the
%                                 opposite trajectory
%
% Laurent Ntibarikure
function [trajTheta, trajPhi, oppTrajTheta, oppTrajPhi] = ...
  getSpiralingHelicoidalTrajectory(N, M, plot, fact)

figProp = getFigureProperties();
map = getColorMap();
tEll=linspace(0,.5,N);
phiEll=tEll*16*pi; 
if nargin>3
  thetaEll=tEll*pi*fact;
else
  thetaEll=tEll*pi;
end
xEll=cos(phiEll).*sin(thetaEll);
yEll=sin(phiEll).*sin(thetaEll); 
zEll=cos(thetaEll);

tEllinv=linspace(0,.5,M);
phiEllinv=tEllinv*16*pi+pi; thetaEllinv=tEllinv*pi;
xEllinv=cos(phiEllinv).*sin(thetaEllinv);
yEllinv=sin(phiEllinv).*sin(thetaEllinv); 
zEllinv=cos(thetaEllinv);
if plot
  figure; plot3(xEll,yEll,zEll, '-', 'LineWidth', 2,...figProp.lw, ...
    'Color', map);
  if ~M
    xlabel('x', 'FontSize', figProp.fs);
    ylabel('y', 'FontSize', figProp.fs);
    zlabel('z', 'FontSize', figProp.fs);
%     legend('Selection trajectory');
  else
    hold on; plot3(xEllinv,yEllinv,zEllinv, '.r', 'LineWidth', figProp.lw);
    xlabel('x', 'FontSize', figProp.fs);
    ylabel('y', 'FontSize', figProp.fs);
    zlabel('z', 'FontSize', figProp.fs);
    legend('Selection trajectory', 'Tested trajectory');
  end
  axis('equal');
  axis([-1 1 -1 1 0 1]);
end
trajTheta = rad2deg(thetaEll);
trajPhi = mod(rad2deg(phiEll),360);
oppTrajTheta = rad2deg(thetaEllinv);
oppTrajPhi = mod(rad2deg(phiEllinv),360);