function plotPolarFF(thetaFF, span, dSpan, pol, gain, gaint, gainp)
    
if pol==1
    maxGain=max(max(gain));
    rho1 = (10*log10(gain(1,:)/maxGain) + span);
    rho1(rho1(:,:)<0)=0;
    rho2 = (10*log10(gain(2,:)/maxGain) + span);
    rho2(rho2(:,:)<0)=0;
    rmax=max(max(rho1), max(rho2));
    figure;
    diagram(rmax,span,dSpan);
    x1 = rho1.*cos(pi/2-thetaFF(1,:));
    y1 = rho1.*sin(pi/2-thetaFF(1,:));
    x2 = rho2.*cos(pi/2-thetaFF(2,:));
    y2 = rho2.*sin(pi/2-thetaFF(2,:));
    plot(x1,y1,'b',x2,y2,'-.r','LineWidth',1.5);
    legend('\phi_{FF}=0°','\phi_{FF}=90°','Location','SouthEast');
    title({['Max. Gain = ',num2str(10*log10(maxGain)),' [dBi]']});
    axis('equal'); axis('off');
elseif pol==2
    maxGain=max(max(gain));
    maxGaint=max(max(gaint));
    maxGainp=max(max(gainp));
    rho1 = (10*log10(gaint(1,:)/maxGain) + span);
    rho1(rho1(:,:)<0)=0;
    rho2 = (10*log10(gaint(2,:)/maxGain) + span);
    rho2(rho2(:,:)<0)=0;

    rho3 = (10*log10(gainp(1,:)/maxGain) + span);
    rho3(rho3(:,:)<0)=0;
    rho4 = (10*log10(gainp(2,:)/maxGain) + span);
    rho4(rho4(:,:)<0)=0;
    rmax=max(max(max(rho1), max(rho2)),max(max(rho3), max(rho4)));
    x1 = rho1.*cos(pi/2-thetaFF(1,:));
    y1 = rho1.*sin(pi/2-thetaFF(1,:));
    x2 = rho2.*cos(pi/2-thetaFF(2,:));
    y2 = rho2.*sin(pi/2-thetaFF(2,:));
    x3 = rho3.*cos(pi/2-thetaFF(1,:));
    y3 = rho3.*sin(pi/2-thetaFF(1,:));
    x4 = rho4.*cos(pi/2-thetaFF(2,:));
    y4 = rho4.*sin(pi/2-thetaFF(2,:));
    figure;
    subplot(1,2,1);
    diagram(rmax,span,dSpan);
    hold on;
    plot(x1,y1,'b',x2,y2,'-.r','LineWidth',1.5);
    legend('\phi_{FF}=0°','\phi_{FF}=90°','Location','SouthEast');
    title({['Max. Gain \theta_{pol} = ',num2str(10*log10(maxGaint)),...
        ' [dBi]']});
    axis('equal'); axis('off');
%         figure;
    subplot(1,2,2);
    diagram(rmax,span,dSpan);
    hold on;
    plot(x3,y3,'b',x4,y4,'-.r','LineWidth',1.5);
    legend('\phi_{FF}=0°','\phi_{FF}=90°','Location','SouthEast');
    title({['Max. Gain \phi_{pol} = ',num2str(10*log10(maxGainp)),...
        ' [dBi]']});
    axis('equal'); axis('off');
end

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
line(rmax*cs,rmax*sn,'linestyle',':','color','k','linewidth',1,...
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
unCenter = -span/(dSpan);
x=(span-3)*cos(pi/2-th1); y=(span-3)*sin(pi/2-th1);
plot(x,y,':k','LineWidth',1,'handlevisibility','off');

text(x(1),y(1)+unCenter,'-3 dB','verticalalignment','bottom','FontSize',8);
while i*dSpan<rmax
    x=(span-i*dSpan)*cos(pi/2-th1); y=(span-i*dSpan)*sin(pi/2-th1);
    plot(x,y,':k','LineWidth',.1,'handlevisibility','off');
    text(x(1),y(1)+unCenter,num2str(-i*dSpan),'verticalalignment','bottom',...
        'FontSize',8);
    i=i+1;
end