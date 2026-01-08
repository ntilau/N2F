% Plots the radiation pattern for XZ and YZ cut-planes in polar diagram
%
% handle = vf_plotFFPolarCutPlanes(handle, span, dSpan, ...
%   thetaFF, gaint, gainp, legend1, ...
%   tRefOrTitle, refGaint, refGainp, legend2, argTitle)
%
% IN: handle = vector figure handles of plot updates (max 2 handles)
%     span = pattern amplitude range in [dB](pattern normalized to the 
%            maximum gain)
%     dSpan = sub-interval range
%     thetaFF = theta look angles for pattern plot
%     gaint-gainp = gain related to theta and phi polarizations
%     legend1 = legend of the pattern
%     tRefOrTitle = title of the diagram or, if OptIN are used, theta for
%                   the refence pattern
%     refGaint-refGainp = refence pattern
%     legend2 = legend of the reference pattern
%     argTitle = title for both patterns
%
% Laurent Ntibarikure
function handle = vf_plotFFPolarCutPlanes(handle, span, dSpan, ...
  thetaFF, gaint, gainp, legend1, ...
  tRefOrTitle, refGaint, refGainp, legend2, argTitle)

  ext=.1; % extention of the bounding circle
% ----- check arguments
  if size(handle,2)<3
    pol = size(handle,2);
  else
    pol = 1;
    disp('Too many handles !!!');
    for i=2:size(handle,2)
      close(handle(i));
    end
  end
  
  plotRef = false;
  insertTitle = false;
  if nargin < 7
    error('Plotting infos incomplete!!!');
  else
    if ~ischar(legend1)
      error('Forgot the legend!!!');
    end
  end
  if nargin == 8
    if ischar(tRefOrTitle)
      insertTitle = true;
      title = tRefOrTitle;
    end
  end
  if nargin > 9
    if nargin < 11
      error('Plotting infos incomplete!!!');
    end
    if ~ischar(legend1) || ~ischar(legend2)
      error('Forgot the legend!!!');
    end
    if nargin > 11
      if ischar(argTitle)
        insertTitle = true;
        title = argTitle;
      end
    end
    plotRef = true;
  end
  
% ----- plots
  figure(handle(1));
  clf;
  gain=gaint+gainp;

  if pol==1 % total fields pattern
    
    if plotRef
      refGain=refGaint+refGainp;
      maxRefGain=max(max(refGain));
      rho3 = (10*log10(refGain(:,1)/maxRefGain) + span);
      rho3(rho3(:,:)<0)=0;
      rho4 = (10*log10(refGain(:,2)/maxRefGain) + span);
      rho4(rho4(:,:)<0)=0;
      rmaxRef=max(max(rho3), max(rho4))+ext;
      x3 = rho3.*cos(pi/2-tRefOrTitle(:,1));
      y3 = rho3.*sin(pi/2-tRefOrTitle(:,1));
      x4 = rho4.*cos(pi/2-tRefOrTitle(:,2));
      y4 = rho4.*sin(pi/2-tRefOrTitle(:,2));
    end
    maxGain=max(max(gain));
    rho1 = (10*log10(gain(:,1)/maxGain) + span);
    rho1(rho1(:,:)<0)=0;
    rho2 = (10*log10(gain(:,2)/maxGain) + span);
    rho2(rho2(:,:)<0)=0;
    rmax=max(max(rho1), max(rho2))+ext;
    if plotRef
      rmax=max(rmax,rmaxRef);
    end
    diagram(rmax,span,dSpan);
    x1 = rho1.*cos(pi/2-thetaFF(:,1));
    y1 = rho1.*sin(pi/2-thetaFF(:,1));
    x2 = rho2.*cos(pi/2-thetaFF(:,2));
    y2 = rho2.*sin(pi/2-thetaFF(:,2));
    if plotRef
      plot(x1,y1,'b',x2,y2,'-.r','LineWidth',1.5);
      hold on;
      plot(x3,y3,'ob',x4,y4,'sr', 'MarkerSize',5);
      legend([legend1,' \phi=0°'],[legend1,' \phi=90°'],...
        [legend2,' \phi=0°'],[legend2,' \phi=90°'],...
        'Location','SouthEast');
    else
      plot(x1,y1,'b',x2,y2,'-.r','LineWidth',1.5);
      legend([legend1,' \phi=0°'],[legend1,' \phi=90°'],...
        'Location','SouthEast');
    end
    text(-1.2*rmax,1.1*rmax,['Max ', ...
      sprintf('%3.4g',10*log10(maxGain)),' [dBi]']);
    axis('equal'); axis('off');
    if insertTitle
      text(.5*rmax,1.1*rmax, title);
    end      
    
  elseif pol==2 % split polarizations
    
    maxGain=max(max(gain));
    maxGaint=max(max(gaint));
    maxGainp=max(max(gainp));
    if plotRef
      refGain=refGaint+refGainp;
      maxRefGain=max(max(refGain));
      rho5 = (10*log10(refGaint(:,1)/maxRefGain) + span);
      rho5(rho5(:,:)<0)=0;
      rho6 = (10*log10(refGaint(:,2)/maxRefGain) + span);
      rho6(rho6(:,:)<0)=0;
      rmaxReft=max(max(rho5), max(rho6))+ext;
      x5 = rho5.*cos(pi/2-tRefOrTitle(:,1));
      y5 = rho5.*sin(pi/2-tRefOrTitle(:,1));
      x6 = rho6.*cos(pi/2-tRefOrTitle(:,2));
      y6 = rho6.*sin(pi/2-tRefOrTitle(:,2));
      rho7 = (10*log10(refGainp(:,1)/maxRefGain) + span);
      rho7(rho7(:,:)<0)=0;
      rho8 = (10*log10(refGainp(:,2)/maxRefGain) + span);
      rho8(rho8(:,:)<0)=0;
      rmaxRefp=max(max(rho7), max(rho8))+ext;
      x7 = rho7.*cos(pi/2-tRefOrTitle(:,1));
      y7 = rho7.*sin(pi/2-tRefOrTitle(:,1));
      x8 = rho8.*cos(pi/2-tRefOrTitle(:,2));
      y8 = rho8.*sin(pi/2-tRefOrTitle(:,2));
    end
    rho1 = (10*log10(gaint(:,1)/maxGain) + span);
    rho1(rho1(:,:)<0)=0;
    rho2 = (10*log10(gaint(:,2)/maxGain) + span);
    rho2(rho2(:,:)<0)=0;
    rho3 = (10*log10(gainp(:,1)/maxGain) + span);
    rho3(rho3(:,:)<0)=0;
    rho4 = (10*log10(gainp(:,2)/maxGain) + span);
    rho4(rho4(:,:)<0)=0;
    rmax=max(max(max(rho1), max(rho2)),max(max(rho3), max(rho4)))+ext;
    x1 = rho1.*cos(pi/2-thetaFF(:,1));
    y1 = rho1.*sin(pi/2-thetaFF(:,1));
    x2 = rho2.*cos(pi/2-thetaFF(:,2));
    y2 = rho2.*sin(pi/2-thetaFF(:,2));
    x3 = rho3.*cos(pi/2-thetaFF(:,1));
    y3 = rho3.*sin(pi/2-thetaFF(:,1));
    x4 = rho4.*cos(pi/2-thetaFF(:,2));
    y4 = rho4.*sin(pi/2-thetaFF(:,2));
    if plotRef
      rmax=max(rmax,rmaxReft);
    end
    diagram(rmax,span,dSpan);
    hold on;
    if plotRef
      plot(x1,y1,'b',x2,y2,'-.r','LineWidth',1.5);
      hold on;
      plot(x5,y5,'ob',x6,y6,'sr', 'MarkerSize',5);
      legend([legend1,' \phi=0°'],[legend1,' \phi=90°'],...
        [legend2,' \phi=0°'],[legend2,' \phi=90°'],...
        'Location','SouthEast');
    else
      plot(x1,y1,'b',x2,y2,'-.r','LineWidth',1.5);
      legend([legend1,' \phi=0°'],[legend1,' \phi=90°'],...
        'Location','SouthEast');
    end
    text(-1.2*rmax,1.1*rmax, ['Max \theta_{pol} ', ...
      sprintf('%3.4g',10*log10(maxGaint)), ' [dBi]']);
    axis('equal'); axis('off');
    if insertTitle
      text(.5*rmax,1.1*rmax, title);
    end  
    
    figure(handle(2));
    clf(handle(2));
    if plotRef
      rmax=max(rmax,rmaxRefp);
    end
    diagram(rmax,span,dSpan);
    hold on;
    if plotRef
      plot(x3,y3,'b',x4,y4,'-.r','LineWidth',1.5);
      hold on;
      plot(x7,y7,'ob',x8,y8,'sr', 'MarkerSize',5);
      legend([legend1,' \phi=0°'],[legend1,' \phi=90°'],...
        [legend2,' \phi=0°'],[legend2,' \phi=90°'],...
        'Location','SouthEast');
    else
      plot(x3,y3,'b',x4,y4,'-.r','LineWidth',1.5);
      legend([legend1,' \phi=0°'],[legend1,' \phi=90°'],...
        'Location','SouthEast');
    end
    text(-1.2*rmax,1.1*rmax,['Max \phi_{pol} ', ...
      sprintf('%3.4g',10*log10(maxGainp)), ' [dBi]']);
    axis('equal'); axis('off');
    if insertTitle
      text(.5*rmax,1.1*rmax, title);
    end
    
  end
  
end

% builds the polar diagram
function diagram(rmax,span,dSpan)
  th1 = linspace(0,2*pi,101);
  xunit = cos(th1);
  yunit = sin(th1);
  inds = 1:(length(th1)-1)/4:length(th1);
  xunit(inds(2:2:4)) = zeros(2,1);
  yunit(inds(1:2:5)) = zeros(3,1);
  patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
    'edgecolor',[0 0 0],'facecolor',[1 1 1],...
    'facealpha', 1,'LineWidth',1, 'handlevisibility','off');
  hold on;
  th2 = (1:6)*2*pi/12;
  cst = cos(th2); snt = sin(th2);
  cs = [-cst; cst];
  sn = [-snt; snt];
  line(rmax*cs,rmax*sn,'linestyle','-','color','k','linewidth',.1,...
    'handlevisibility','off');
  rt = 1.1*rmax;
  for i = 1:length(th2)
    text(rt*cst(i),rt*snt(i),[int2str((90-i*30)),...
      '°'], 'horizontalalignment','center',...
      'handlevisibility','off');
    if i == length(th2)
      loc = int2str(90);
    elseif i ~= 3
      val = - 90- sign(i-6)*i*30;
      loc = int2str( -(- sign(val)*180 + val));
    else
      loc = int2str( 180 );
    end
    text(-rt*cst(i),-rt*snt(i),[loc '°'],'horizontalalignment','center',...
    'handlevisibility','off')
  end
  hold on;
  i=1;
  unCenter = 0;
  x=(span-3)*cos(pi/2-th1); y=(span-3)*sin(pi/2-th1);
  plot(x,y,'-k','LineWidth',.1,'handlevisibility','off');
  text(x(1)+.5,y(1)+unCenter,'-3 dB','verticalalignment','bottom','FontSize',8);
  while i*dSpan<rmax
    x=(span-i*dSpan)*cos(pi/2-th1); y=(span-i*dSpan)*sin(pi/2-th1);
    plot(x,y,'-k','LineWidth',.1,'handlevisibility','off');
    text(x(1)+.5,y(1)+unCenter,num2str(-i*dSpan),'verticalalignment','bottom',...
      'FontSize',8);
    i=i+1;
  end
end
