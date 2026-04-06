clear; close all; clc;
addpath('../../matlabLib');
%% constants
pi=3.14159265358979323846;
eps0=8.85418781761e-12;
mu0=pi*4e-7;
c0=1/sqrt(eps0*mu0);
zeta0=sqrt(mu0/eps0);
freq=1e9;
lambda0=c0/freq;
k0=2*pi/lambda0;
%% sphere points
rho=1; % [m]
dtheta=2; dphi=2; % [°] high oversampling to get a good field's accuracy
theta=(dtheta/2:dtheta:180-dtheta/2)*pi/180;
phi=(0:dphi:360-dphi)*pi/180;
[theta,phi]=meshgrid(theta,phi);
[rhox,rhoy,rhoz]=sph2cart(phi,pi/2-theta,rho);
[rx,ry,rz]=sph2cart(phi,pi/2-theta,1);
dS=rho.^2.*sin(theta).*(2*pi^2/(size(theta,2)*size(phi,1)));
%% electric/magnetic dipole location and orientation (theta, phi)
xd=0; yd=0; zd=0; % position [m]
tJ=90; pJ=270; tM=90; pM=0; % direction (theta,phi)
[Jx,Jy,Jz]=spherical2cartesian(1,0,0,tJ*pi/180,pJ*pi/180); % 1 [Am] electric dipole
[Mx,My,Mz]=spherical2cartesian(zeta0,0,0,tM*pi/180,pM*pi/180); % 1 [Vm] magnetic dipole
%% R
Rx=rhox-xd; Ry=rhoy-yd; Rz=rhoz-zd;
R=sqrt(Rx.^2+Ry.^2+Rz.^2);
RxV=Rx./R; RyV=Ry./R; RzV=Rz./R;
%% near field computations from Milligan p.46
JdotRV=Jx.*RxV+Jy.*RyV+Jz.*RzV;
[JcrossRVx,JcrossRVy,JcrossRVz]=crossOperator(Jx,Jy,Jz,RxV,RyV,RzV);
MdotRV=Mx.*RxV+My.*RyV+Mz.*RzV;
[McrossRVx,McrossRVy,McrossRVz]=crossOperator(Mx,My,Mz,RxV,RyV,RzV);
Ex=( zeta0*k0^2/(4*pi).*( Jx.*(-1i./(k0*R)-1./(k0*R).^2+1i./(k0*R).^3) +...
    JdotRV.*RxV.*(1i./(k0*R)+3./(k0*R).^2-3i./(k0*R).^3) ) - k0^2/(4*pi).*...
    McrossRVx.*(1i./(k0.*R)+1./(k0.*R).^2) ).*exp(-1i*k0.*R);
Ey=( zeta0*k0^2/(4*pi).*( Jy.*(-1i./(k0*R)-1./(k0*R).^2+1i./(k0*R).^3) +...
    JdotRV.*RyV.*(1i./(k0*R)+3./(k0*R).^2-3i./(k0*R).^3) ) - k0^2/(4*pi).*...
    McrossRVy.*(1i./(k0.*R)+1./(k0.*R).^2) ).*exp(-1i*k0.*R);
Ez=( zeta0*k0^2/(4*pi).*( Jz.*(-1i./(k0*R)-1./(k0*R).^2+1i./(k0*R).^3) +...
    JdotRV.*RzV.*(1i./(k0*R)+3./(k0*R).^2-3i./(k0*R).^3) ) - k0^2/(4*pi).*...
    McrossRVz.*(1i./(k0.*R)+1./(k0.*R).^2) ).*exp(-1i*k0.*R);
Hx=( k0^2/(4*pi*zeta0).*( Mx.*(-1i./(k0*R)-1./(k0*R).^2+1i./(k0*R).^3) +...
    MdotRV.*RxV.*(1i./(k0*R)+3./(k0*R).^2-3i./(k0*R).^3) ) + k0^2/(4*pi).*...
    JcrossRVx.*(1i./(k0.*R)+1./(k0.*R).^2) ) .*exp(-1i*k0.*R);
Hy=( k0^2/(4*pi*zeta0).*( My.*(-1i./(k0*R)-1./(k0*R).^2+1i./(k0*R).^3) +...
    MdotRV.*RyV.*(1i./(k0*R)+3./(k0*R).^2-3i./(k0*R).^3) ) + k0^2/(4*pi).*...
    JcrossRVy.*(1i./(k0.*R)+1./(k0.*R).^2) ).*exp(-1i*k0.*R);
Hz=( k0^2/(4*pi*zeta0).*( Mz.*(-1i./(k0*R)-1./(k0*R).^2+1i./(k0*R).^3) +...
    MdotRV.*RzV.*(1i./(k0*R)+3./(k0*R).^2-3i./(k0*R).^3) ) + k0^2/(4*pi).* ...
    JcrossRVz.*(1i./(k0.*R)+1./(k0.*R).^2) ).*exp(-1i*k0.*R);
[Er,Et,Ep]=cartesian2spherical(Ex,Ey,Ez,theta,phi);
[Hr,Ht,Hp]=cartesian2spherical(Hx,Hy,Hz,theta,phi);
%% field plots
if(0)
    extrusionParam = sqrt(abs(Er).^2+abs(Et).^2+abs(Ep).^2);
    [x,y,z] = sph2cart(phi,pi/2-theta,extrusionParam);
    figure; surf(x,y,z,extrusionParam,'EdgeColor','none');
    axis('equal');%view(0,0);
end
%% equivalent currents
[Jsx,Jsy,Jsz] = crossOperator(-rx,-ry,-rz,Hx,Hy,Hz);
[Msx,Msy,Msz] = crossOperator(rx,ry,rz,Ex,Ey,Ez);
%% nf total power radiated
% [Sr,St,Sp]=crossOperator(Er,Et,Ep,conj(Hr),conj(Ht),conj(Hp));
% Pr=1/2*real(sum(sum(Sr.*dS)));
% fprintf('Radiated power : %g Watts\n',Pr);
[Sx,Sy,Sz]=crossOperator(Jsx,Jsy,Jsz,conj(Msx),conj(Msy),conj(Msz));
Pr=1/2*real(sum(sum((Sx.*rx+Sy.*ry+Sz.*rz).*dS)));
fprintf('Radiated power : %g Watts\n',Pr);
%%
if(0)
  figure;
  surf(rhox,rhoy,rhoz,sqrt(abs(Ex).^2+abs(Ey).^2+abs(Ez).^2),...
      'FaceAlpha',.5,'EdgeAlpha',.5,'EdgeColor','none'); hold on;
  quiver3(rhox,rhoy,rhoz,real(Jsx),real(Jsy),real(Jsz),1,'r');
  quiver3(rhox,rhoy,rhoz,real(Msx),real(Msy),real(Msz),1,'b');
  axis('equal'); view(0,0);
end
%% NF2FF
ffPlot='3D';
dthetaFF=5; dphiFF=5;
switch ffPlot
case '3D'
  thetaFF=(0:dthetaFF:180)*pi/180;
  phiFF=(0:dphiFF:360)*pi/180;
  [thetaFF,phiFF]=meshgrid(thetaFF,phiFF);
  ExFF=zeros(size(thetaFF));
  EyFF=ExFF;EzFF=ExFF;
  tic
  for i=1:size(thetaFF,2)
    for j=1:size(phiFF,1)
%       fprintf('thetaFF=%g, phiFF=%g\n',thetaFF(j,i)*180/pi,phiFF(j,i)*180/pi);
      RffVx=sin(thetaFF(j,i)).*cos(phiFF(j,i));
      RffVy=sin(thetaFF(j,i)).*sin(phiFF(j,i));
      RffVz=cos(thetaFF(j,i));
      green=exp(1i*k0*(RffVx.*rhox+RffVy.*rhoy+RffVz.*rhoz));
      JsdotRffV=Jsx.*RffVx+Jsy.*RffVy+Jsz.*RffVz;
      [MscrossRffVx,MscrossRffVy,MscrossRffVz]=...
          crossOperator(Msx,Msy,Msz,RffVx,RffVy,RffVz);
      ExFF(j,i)=-1i*k0/(4*pi).*sum(sum( ...
          (zeta0*(Jsx-JsdotRffV.*RffVx) + MscrossRffVx).*green.*dS ));
      EyFF(j,i)=-1i*k0/(4*pi).*sum(sum( ...
          (zeta0*(Jsy-JsdotRffV.*RffVy) + MscrossRffVy).*green.*dS ));
      EzFF(j,i)=-1i*k0/(4*pi).*sum(sum( ...
          (zeta0*(Jsz-JsdotRffV.*RffVz) + MscrossRffVz).*green.*dS ));
    end
  end
  fprintf('nf2ff computation time : %g s.\n',toc)
case '2D'
  thetaFF=(-180:1:180)*pi/180;
  phiFF=[0 90]*pi/180;
  [thetaFF,phiFF]=meshgrid(thetaFF,phiFF);
  ExFF=zeros(size(thetaFF));
  EyFF=ExFF;EzFF=ExFF;
  tic
  for i=1:size(thetaFF,2)
    for j=1:size(phiFF,1)
      RffVx=sin(thetaFF(j,i)).*cos(phiFF(j,i));
      RffVy=sin(thetaFF(j,i)).*sin(phiFF(j,i));
      RffVz=cos(thetaFF(j,i));
      green=exp(1i*k0*(RffVx.*rhox+RffVy.*rhoy+RffVz.*rhoz));
      JsdotRffV=Jsx.*RffVx+Jsy.*RffVy+Jsz.*RffVz;
      [MscrossRffVx,MscrossRffVy,MscrossRffVz]=...
          crossOperator(Msx,Msy,Msz,RffVx,RffVy,RffVz);
      ExFF(j,i)=-1i*k0/(4*pi).*sum(sum( ...
          (zeta0*(Jsx-JsdotRffV.*RffVx) + MscrossRffVx).*green.*dS ));
      EyFF(j,i)=-1i*k0/(4*pi).*sum(sum( ...
          (zeta0*(Jsy-JsdotRffV.*RffVy) + MscrossRffVy).*green.*dS ));
      EzFF(j,i)=-1i*k0/(4*pi).*sum(sum( ...
          (zeta0*(Jsz-JsdotRffV.*RffVz) + MscrossRffVz).*green.*dS ));
    end
  end
  fprintf('nf2ff computation time : %g s.\n',toc)
end
%% electric FF spherical components in [V] (exp(-jkR)/R term neglected)
[ErFF,EtFF,EpFF]=cartesian2spherical(ExFF,EyFF,EzFF,thetaFF,phiFF);
%% Directivity or ideal Gain
reSrFF=1/2*(abs(EtFF).^2+abs(EpFF).^2)./zeta0; % real part of Poynting vector
Gain=4*pi/Pr.*reSrFF; % Directivity
GainVP=2*pi/Pr.*abs(EtFF).^2./zeta0; % Directivity vertical polarization
GainHP=2*pi/Pr.*abs(EpFF).^2./zeta0; % Directivity horizontal polarization
%% radiation solid
Polarization='Split';
if strcmp(ffPlot,'3D')
  switch Polarization
    case 'Total'
      [x,y,z] = sph2cart(phiFF,pi/2-thetaFF,Gain);
      figure; surf(x,y,z,Gain,'FaceAlpha',.75,'EdgeAlpha',.5,...
          'EdgeColor','none'); hold on;
      quiver3(x,y,z,real(ExFF),real(EyFF),real(EzFF),1,'k');
      axis('equal');xlabel('x');ylabel('y');zlabel('z');%view(0,0);
      title('Radiation solid - Total fields');
    case 'Split'
      figure;
      % vertical
      subplot(1,2,1);
      [x,y,z] = sph2cart(phiFF,pi/2-thetaFF,GainVP);
      surf(x,y,z,GainVP,'FaceAlpha',.75,'EdgeAlpha',.5,...
          'EdgeColor','none'); hold on;
      quiver3(x,y,z,real(ExFF),real(EyFF),real(EzFF),1,'k');
      axis('equal');xlabel('x');ylabel('y');zlabel('z');%view(0,0);
      title('\theta polarization');
      % horizontal
      subplot(1,2,2);
      [x,y,z] = sph2cart(phiFF,pi/2-thetaFF,GainHP);
      surf(x,y,z,GainHP,'FaceAlpha',.75,'EdgeAlpha',.5,...
          'EdgeColor','none'); hold on;
      quiver3(x,y,z,real(ExFF),real(EyFF),real(EzFF),1,'k');
      axis('equal');xlabel('x');ylabel('y');zlabel('z');%view(0,0);
      title('\phi polarization');
  end
else
  maxGain=max(max(Gain));
  figure;
  plot(thetaFF(1,:)*180/pi, 10*log10(Gain(1,:)),'b',...
      thetaFF(2,:)*180/pi, 10*log10(Gain(2,:)),'-.r');
  axis('tight');v=axis;axis([v(1) v(2) -60 v(4)]);
  xlabel('\theta_{FF}');ylabel('Gain [dB]');
  legend('\phi_{FF}=0°','\phi_{FF}=90°','Location','SouthEast');
  title(['Maximum Gain = ', num2str(10*log10(maxGain)),' dB']);
end
%% used to compare the results with FEKO simulation (HertzDipole.txt)
checkMagEt=abs(EtFF);
checkAngEt=angle(EtFF)*180/pi;
checkMagEp=abs(EpFF);
checkAngEp=angle(EpFF)*180/pi;
checkGain=10*log10(Gain);