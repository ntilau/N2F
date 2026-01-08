function plotFF(thetaFF, span, pol, gaint, gainp)

  gain = gaint+gainp;
  maxGain=max(max(gain));
  if pol==1
      scf();
      plot(thetaFF(1,:)*180/%pi, 10*log10(gain(1,:)),'b',...
      thetaFF(2,:)*180/%pi, 10*log10(gain(2,:)),'-.r');
      a=gca();
      a.tight_limits='on';
      a.data_bounds(1:2,2)= [-span; 10*log10(maxGain)];
      xtitle(['Max. gain = ', sprintf('%2.4g',10*log10(maxGain)) ,' [dBi]'],...
          'theta','gain [dBi]');
      legend('phi=0','phi=90',4);
  elseif pol==2
      maxGaint=max(max(gaint));
      maxGainp=max(max(gainp));
      scf();
      plot(thetaFF(1,:)*180/%pi, 10*log10(gaint(1,:)),'b',...
      thetaFF(2,:)*180/%pi, 10*log10(gaint(2,:)),'-.r');
      a=gca();
      a.tight_limits='on';
      a.data_bounds(1:2,2)= [-span; 10*log10(maxGain)];
      xlabel('theta');ylabel('gain [dBi]');
      legend('phi=0','phi=90',4);
      title(['Max. gain theta_pol = ', sprintf('%2.4g',10*log10(maxGaint)),'[dBi]']);
      
      scf();
      plot(thetaFF(1,:)*180/%pi, 10*log10(gainp(1,:)),'b',...
          thetaFF(2,:)*180/%pi, 10*log10(gainp(2,:)),'-.r');
      a2=gca();
      a2.tight_limits='on';
      a2.data_bounds(1:2,2)= [-span; 10*log10(maxGain)];
      xlabel('theta');ylabel('gain [dBi]');
      legend('phi=0','phi=90',4);
      title(['Max. gain phi_pol = ', sprintf('%2.4g',10*log10(maxGainp)),' [dBi]']);
  end

endfunction