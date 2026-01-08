//% Computes the far fields from surface currents on the box
function [gaint, gainp, thetaFF, phiFF] =  ...
  farFieldSolver(k0, z0, boxPos, dS, dthetaFF, Js, Ms, Pr)

  exec('cartesian2spherical.sci',-1);
  exec('crossOperator.sci',-1);

  thetaFF=(-180:dthetaFF:180)*%pi/180;
  phiFF=[0 90]*%pi/180;
  [thetaFF,phiFF]=meshgrid(thetaFF,phiFF);
  ExFF=zeros(thetaFF);
  EyFF=ExFF;EzFF=ExFF;
  tic();
  for i=1:size(thetaFF,2)
      for j=1:size(phiFF,1)
          RffVx=sin(thetaFF(j,i)).*cos(phiFF(j,i));
          RffVy=sin(thetaFF(j,i)).*sin(phiFF(j,i));
          RffVz=cos(thetaFF(j,i));
  
          green=exp(%i*k0*(RffVx.*boxPos(1,:)+...
              RffVy.*boxPos(2,:)+RffVz.*boxPos(3,:)));

          JsdotRffV=Js(1,:).*RffVx+...
              Js(2,:).*RffVy+Js(3,:).*RffVz;
          [MscrossRffVx,MscrossRffVy,MscrossRffVz]=...
              crossOperator(Ms(1,:),...
              Ms(2,:),Ms(3,:),RffVx,RffVy,RffVz);
          ExFF(j,i)= ExFF(j,i) -%i*k0/(4*%pi).*sum(sum( ...
              (z0*(Js(1,:)-JsdotRffV.*RffVx) +...
              MscrossRffVx).*green.*dS ));
          EyFF(j,i)= EyFF(j,i) -%i*k0/(4*%pi).*sum(sum( ...
              (z0*(Js(2,:)-JsdotRffV.*RffVy) +...
              MscrossRffVy).*green.*dS ));
          EzFF(j,i)= EzFF(j,i) -%i*k0/(4*%pi).*sum(sum( ...
              (z0*(Js(3,:)-JsdotRffV.*RffVz) +...
              MscrossRffVz).*green.*dS ));
  
      end
  end
  printf('##> nf2ff computation time : %g s.\n',toc())
  // electric FF spherical components in [V] (exp(-jkR)/R term neglected)
  [ErFF,EtFF,EpFF]=cartesian2spherical(ExFF,EyFF,EzFF,thetaFF,phiFF);
  // Directivity or ideal Gain
  gaint=2*%pi/Pr.*abs(EtFF).^2./z0;
  gainp=2*%pi/Pr.*abs(EpFF).^2./z0;

endfunction