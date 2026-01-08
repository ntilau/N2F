% Plots the far field pattern on the cut planes (phi=cost) selected
%
% vf_plotFFCutPlanes(thetaFF, gaint, gainp, span, pol,...
%   thetaFFref, refGaint, refGainp)
%
% IN: thetaFF = theta look angles on the constant phi plane
%     gaint = gain related to theta polarized electric field in [dBi]
%     gainp = gain related to phi polarized electric field in [dBi]
%     span = pattern range in [dB]
%     pol = polarizations 1: total, 2: split in theta/phi
%     thetaFFref = theta look angles on the constant phi plane for the
%                reference values - allows decimation
%     refGaint = gain related to theta polarized electric field in [dBi] of 
%                the reference
%     refGainp = gain related to phi polarized electric field in [dBi] of 
%                the reference
% OptIN: matchedT = vertical lines for selected angles in [rad]
%
% Laurent Ntibarikure
function vf_plotFFCutPlanes(thetaFF, gaint, gainp, span, pol,...
  thetaFFref, refGaint, refGainp, matchedT)

figProp.ms = 3;
figProp.lw = 1.5;
figProp.fs = 12;

plotRef = true;
if nargin < 6
  plotRef = false;
end

gain=gaint+gainp;
maxGain=max(max(gain));
if plotRef
  refGain=refGaint+refGainp;
  maxRefGain=max(max(refGain));
end
a(1)=min(min(thetaFF))*180/pi;
a(2)=max(max(thetaFF))*180/pi;

switch pol
  
  case 1
  figure;
  plot(rad2deg(thetaFF(:,1)), 10*log10(gain(:,1)),'b',...
    rad2deg(thetaFF(:,2)), 10*log10(gain(:,2)),'-.r', ...
    'MarkerSize', figProp.ms,'LineWidth', figProp.lw);
  if plotRef
    hold on;
    plot(rad2deg(thetaFFref(:,1)), 10*log10(refGain(:,1)),'*b',...
      rad2deg(thetaFFref(:,2)), 10*log10(refGain(:,2)),'+r',...
      'MarkerSize', figProp.ms,'LineWidth', figProp.lw);
    a(3) = max(10*log10(maxRefGain),10*log10(maxGain));
    axis([a(1) a(2) -span a(3)]);
    if nargin >8
      for i=1:length(matchedT)
        line([matchedT(i) matchedT(i)]*180/pi, [-span a(3)], ...
          'LineStyle', '-.', 'Color','k', 'LineWidth', figProp.lw);
      end
    end
    xlabel('\theta_{FF} [°]', 'FontSize', figProp.fs);
    ylabel('gain [dBi]', 'FontSize', figProp.fs);
    legend('\phi_{FF}=0°','\phi_{FF}=90°','Location','SouthEast');
    title({['Max ', sprintf('%3.4g',10*log10(maxGain)),' [dBi]'];...
        ['Max ref. ', sprintf('%3.4g',10*log10(maxRefGain)),' [dBi]']}, ...
        'FontSize', figProp.fs);
  else
    a(3) = 10*log10(maxGain);
    axis([a(1) a(2) -span a(3)]);
    if nargin >8
      for i=1:length(matchedT)
        line([matchedT(i) matchedT(i)]*180/pi, [-span a(3)], ...
          'LineStyle', '-.', 'Color','k', 'LineWidth', figProp.lw);
      end
    end
    xlabel('\theta_{FF} [°]', 'FontSize', figProp.fs);     
    ylabel('gain [dBi]', 'FontSize', figProp.fs);
    legend('\phi_{FF}=0°','\phi_{FF}=90°','Location','SouthEast');
    title(['Max ', sprintf('%3.4g',10*log10(maxGain)),' [dBi]'], ...
        'FontSize', figProp.fs);
  end

  case 2
  maxGaint=max(max(gaint));
  maxGainp=max(max(gainp));
  if plotRef
    maxRefGaint=max(max(refGaint));
    maxRefGainp=max(max(refGainp));
  end

  figure;
  plot(rad2deg(thetaFF(:,1)), 10*log10(gaint(:,1)),'b',...
    rad2deg(thetaFF(:,2)), 10*log10(gaint(:,2)),'-.r', ...
    'MarkerSize', figProp.ms,'LineWidth', figProp.lw);
  if plotRef
    hold on;
    plot(rad2deg(thetaFFref(:,1)), 10*log10(refGaint(:,1)),'*b',...
      rad2deg(thetaFFref(:,2)), 10*log10(refGaint(:,2)),'+r',...
      'MarkerSize', figProp.ms,'LineWidth', figProp.lw);
    a(3) = max(10*log10(maxRefGaint),10*log10(maxGaint));
    axis([a(1) a(2) -span a(3)]);
    if nargin >8
      for i=1:length(matchedT)
        line([matchedT(i) matchedT(i)]*180/pi, [-span a(3)], ...
          'LineStyle', '-.', 'Color','k', 'LineWidth', figProp.lw);
      end
    end
    xlabel('\theta_{FF} [°]', 'FontSize', figProp.fs);     
    ylabel('gain [dBi]', 'FontSize', figProp.fs);
    legend('\phi_{FF}=0°','\phi_{FF}=90°','Location','SouthEast');
    title({['Max \theta_{pol} ', ...
      sprintf('%3.4g',10*log10(maxGaint)),' [dBi]'];...
      ['Max ref. \theta_{pol} ', ...
      sprintf('%3.4g',10*log10(maxRefGaint)),' [dBi]']}, ...
        'FontSize', figProp.fs);
  else
    a(3) = 10*log10(maxGaint);
    axis([a(1) a(2) -span a(3)]);
    if nargin >8
      for i=1:length(matchedT)
        line([matchedT(i) matchedT(i)]*180/pi, [-span a(3)], ...
          'LineStyle', '-.', 'Color','k', 'LineWidth', figProp.lw);
      end
    end
    xlabel('\theta_{FF} [°]', 'FontSize', figProp.fs);     
    ylabel('gain [dBi]', 'FontSize', figProp.fs);
    legend('\phi_{FF}=0°','\phi_{FF}=90°','Location','SouthEast');
    title(['Max \phi_{pol} ', ...
      sprintf('%3.4g',10*log10(maxGainp)),' [dBi]'], ...
      'FontSize', figProp.fs);
  end

  figure;
  plot(rad2deg(thetaFF(:,1)), 10*log10(gainp(:,1)),'b',...
    rad2deg(thetaFF(:,2)), 10*log10(gainp(:,2)),'-.r', ...
    'MarkerSize', figProp.ms,'LineWidth', figProp.lw);
  if plotRef
    hold on;
    plot(rad2deg(thetaFFref(:,1)), 10*log10(refGainp(:,1)),'*b',...
      rad2deg(thetaFFref(:,2)), 10*log10(refGainp(:,2)),'+r',...
      'MarkerSize', figProp.ms,'LineWidth', figProp.lw);
    a(3) = max(10*log10(maxRefGainp),10*log10(maxGainp));
    axis([a(1) a(2) -span a(3)]);
    if nargin >8
      for i=1:length(matchedT)
        line([matchedT(i) matchedT(i)]*180/pi, [-span a(3)], ...
          'LineStyle', '-.', 'Color','k', 'LineWidth', figProp.lw);
      end
    end
    xlabel('\theta_{FF} [°]', 'FontSize', figProp.fs);     
    ylabel('gain [dBi]', 'FontSize', figProp.fs);
    legend('\phi_{FF}=0°','\phi_{FF}=90°','Location','SouthEast');
    title({['Max \phi_{pol} ', ...
      sprintf('%3.4g',10*log10(maxGainp)),' [dBi]'];...
      ['Max ref. \phi_{pol} ', ...
      sprintf('%3.4g',10*log10(maxRefGainp)),' [dBi]']}, ...
      'FontSize', figProp.fs);
  else
    a(3) = 10*log10(maxGainp);
    axis([a(1) a(2) -span a(3)]);
    if nargin >8
      for i=1:length(matchedT)
        line([matchedT(i) matchedT(i)]*180/pi, [-span a(3)], ...
          'LineStyle', '-.', 'Color','k', 'LineWidth', figProp.lw);
      end
    end
    xlabel('\theta_{FF} [°]', 'FontSize', figProp.fs);     
    ylabel('gain [dBi]', 'FontSize', figProp.fs);
    legend('\phi_{FF}=0°','\phi_{FF}=90°','Location','SouthEast');
    title(['Max \phi_{pol} ', ...
      sprintf('%3.4g',10*log10(maxGainp)),' [dBi]'], ...
      'FontSize', figProp.fs);
  end
end