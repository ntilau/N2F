% Computes the far electric field detectors from the near fiels computed on 
% a surface using Huygens' principle
%
% [EtFF, EpFF] =  vf_nf2ffSolver(k0, z0, surfPos, n, dS,...
% E, H, thetaFF, phiFF)
%
% IN: k0-z0 = electrical parameters
%     surfPos = near fields sampling points
%     n = normal unit vector outwardly directed
%     dS = surface patches areas
%     E = electric field vector in Cartesian coordinates
%     M = magnetic field vector in Cartesian coordinates
%     thetaFF-phiFF = far field look angles (meshgrid)
%
% OUT: [EtFF, EpFF] = spherical components of the far electic field 
%                     detector for the choosen look angles
%
% Laurent Ntibarikure
function [EtFF, EpFF] =  vf_nf2ffSolver(k0, z0, surfPos, n, dS,...
  E, H, thetaFF, phiFF)

tic();
Js = cross(n,H);
Ms = cross(E, n);
ExFF=zeros(size(thetaFF));
EyFF=ExFF;EzFF=ExFF;
for i=1:size(thetaFF,2)
  for j=1:size(phiFF,1)
    RffVx=sin(thetaFF(j,i)).*cos(phiFF(j,i));
    RffVy=sin(thetaFF(j,i)).*sin(phiFF(j,i));
    RffVz=cos(thetaFF(j,i));
    green=exp(1i*k0*(RffVx.*surfPos(1,:)+ ...
      RffVy.*surfPos(2,:)+RffVz.*surfPos(3,:)));
    JsdotRffV=Js(1,:).*RffVx+ Js(2,:).*RffVy+Js(3,:).*RffVz;
    [MscrossRffVx,MscrossRffVy,MscrossRffVz]=...
      crossOperator(Ms(1,:),Ms(2,:),Ms(3,:),RffVx,RffVy,RffVz);
    % function trapz() does not significantly enhance the integration
    % seeking for surface Gaussian quadrature :-|
    ExFF(j,i)= ExFF(j,i) - 1i*k0/(4*pi).*sum( ...
      (z0*(Js(1,:)-JsdotRffV.*RffVx) + MscrossRffVx).*green.*dS );
    EyFF(j,i)= EyFF(j,i) - 1i*k0/(4*pi).*sum( ...
      (z0*(Js(2,:)-JsdotRffV.*RffVy) + MscrossRffVy).*green.*dS );
    EzFF(j,i)= EzFF(j,i) - 1i*k0/(4*pi).*sum( ...
      (z0*(Js(3,:)-JsdotRffV.*RffVz) + MscrossRffVz).*green.*dS );
  end
end
fprintf('#> nf2ff computation time : %g s.\n',toc)
% electric FF spherical components in [V] (exp(-jkR)/R term neglected)
[ErFF,EtFF,EpFF]=cartesian2spherical(ExFF,EyFF,EzFF,thetaFF,phiFF);