% Plots the radiation solid
%
% vf_plotFF3d(handle, thetaFF, phiFF, offset, gaint, gainp)
%
% IN: handle = figure handle for plot update
%     thetaFF-phiFF = look angles (meshgrid)
%     offset = pattern range in [dB]
%     gaint-gainp = gain for 3d pattern
%     infos = further informations added to the title
%
% Laurent Ntibarikure
function handle = vf_plotFF3d(handle, thetaFF, phiFF, offset, ...
  gaint, gainp, infos)

figure(handle);
clf();
gain=gaint+gainp;
maxGain=max(max(gain));
r=10*log10(gain/maxGain)+offset;
r(r(:,:)<0)=0;
[x,y,z] = sph2cart(phiFF,pi/2-thetaFF,r);
surf(x,y,z,r,'FaceAlpha',1,'EdgeAlpha',.1,...
    'EdgeColor','k'); 
hold on;
axis('equal');
view(135,45);
axis('off'); 
line([0 0],[0 0],[0 3+offset],'color','k');
line([0 0],[0 3+offset],[0 0],'color','k');
line([0 3+offset],[0 0],[0 0],'color','k');
text(0,0,3+offset,'z');
text(0,3+offset,0,'y');
text(3+offset,0,0,'x');
if nargin>6
  title({['Radiation solid. Max ', sprintf('%3.4g',10*log10(maxGain)), ...
    ' [dBi]']; infos});
else
  title({['Radiation solid. Max ', sprintf('%3.4g',10*log10(maxGain)), ...
    ' [dBi]']});
end