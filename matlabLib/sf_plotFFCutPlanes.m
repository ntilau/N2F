% Plots the far field pattern on the cut planes (phi=cost) selected
%
% sf_plotFFCutPlanes(theta, gain, refTheta, refGain, planes,...
%   infos, filename)
%
% IN: theta = theta look angles on the constant phi plane
%     gain = directivity values in [dBi]
%     refTheta = theta look angles on the constant phi plane for the
%                reference values - allows decimation
%     refGain = directivity values in [dBi] of the reference
%     planes = vector of the planes to plot (row indices of gain)
%     infos = vector of structures. For the i-th plot:
%             infos(i).markers = string of chosen marker and color
%             infos(i).title = string of chosen title
%             infos(i).legend1 = string of chosen legend for gain
%             infos(i).legend2 = string of chosen legend for refGain
%     filename{i} = cell of filenames (strings) for printing with printEPS
%
% Laurent Ntibarikure
function sf_plotFFCutPlanes(theta, gain, refTheta, refGain, planes,...
  infos, filename)

figProp = getFigureProperties();

for i=1:length(planes)
  figure;

  if nargin>5
    if isfield(infos(i), 'markers')
      plot(theta(1,:)*180/pi, 10*log10(gain(planes(i),:)), ...
        infos(i).markers{1}, ...
        'LineWidth', figProp.lw, 'MarkerSize', figProp.ms);
      hold on;
      plot(refTheta(1,:)*180/pi, 10*log10(refGain(planes(i),:)), ...
        infos(i).markers{2}, ...
        'LineWidth', figProp.lw, 'MarkerSize', figProp.ms);
    else
      plot(theta(1,:)*180/pi, 10*log10(gain(planes(i),:)), '-b', ...
        'LineWidth', figProp.lw, 'MarkerSize', figProp.ms);
      hold on;
      plot(refTheta(1,:)*180/pi, 10*log10(refGain(planes(i),:)), '+r', ...
        'LineWidth', figProp.lw, 'MarkerSize', figProp.ms);
    end
  else
    plot(theta(1,:)*180/pi, 10*log10(gain(planes(i),:)), '-b', ...
      'LineWidth', figProp.lw, 'MarkerSize', figProp.ms);
    hold on;
    plot(refTheta(1,:)*180/pi, 10*log10(refGain(planes(i),:)), '+r', ...
      'LineWidth', figProp.lw, 'MarkerSize', figProp.ms);
  end
  axis('tight'); v = axis;
  if min(min(refGain))<1e-6, axis([v(1) v(2) -60 v(4)]); end
  line([-90 -90], [v(3) v(4)], 'Color', 'k', 'LineStyle', '-.');
  line([90 90], [v(3) v(4)], 'Color', 'k', 'LineStyle', '-.');
  line([270 270], [v(3) v(4)], 'Color', 'k', 'LineStyle', '-.');
  xlabel('\theta [°]','fontsize', figProp.fs); 
  ylabel('Directivity [dBi]','fontsize', figProp.fs);
  if nargin>5
    if isfield(infos(i), 'title')
      title(infos(i).title);
    end
    if isfield(infos(i), 'legend1')
      if isfield(infos(i), 'legend2')
        legend(infos(i).legend1, infos(i).legend2, 'Location','SouthEast');
      end
    end
  else
    legend('N2F','Direct','Location','SouthEast');
  end
  if nargin>6
    printEPS('',filename{i});
  end
end

