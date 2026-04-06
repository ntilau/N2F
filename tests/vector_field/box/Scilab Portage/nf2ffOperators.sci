// TOO COMPUTATIONALLY EXPENSIVE FOR SCILAB -> USE FarFieldSolver()
//
//% Computes the nf2ff operators for arbitrarily shaped surface
function [Op,thetaFF,phiFF] =  nf2ffOperators(k0, z0, pos, dS, dthetaFF)

  stacksize('max');
  thetaFF=(-180:dthetaFF:180)*%pi/180;
  phiFF=[0 90]*%pi/180;
  [thetaFF,phiFF]=meshgrid(thetaFF,phiFF);
  OpDim = zeros(1,size(pos,2));
  OpxJx=OpDim; OpxJy=OpDim; OpxJz=OpDim;
  OpxMx=OpDim; OpxMy=OpDim; OpxMz=OpDim;
  OpyJx=OpDim; OpyJy=OpDim; OpyJz=OpDim;
  OpyMx=OpDim; OpyMy=OpDim; OpyMz=OpDim;
  OpzJx=OpDim; OpzJy=OpDim; OpzJz=OpDim;
  OpzMx=OpDim; OpzMy=OpDim; OpzMz=OpDim;
  
  OpDim = zeros(size(thetaFF,2),size(pos,2),size(phiFF,1));
  Op = struct('tJx',OpDim,'tJy',OpDim,'tJz',OpDim,...
    'tMx',OpDim,'tMy',OpDim,'tMz',OpDim,...
    'pJx',OpDim,'pJy',OpDim,'pJz',OpDim,...
    'pMx',OpDim,'pMy',OpDim,'pMz',OpDim);

  for i=1:size(thetaFF,2)
      for j=1:size(phiFF,1)
          RffVx=sin(thetaFF(j,i)).*cos(phiFF(j,i));
          RffVy=sin(thetaFF(j,i)).*sin(phiFF(j,i));
          RffVz=cos(thetaFF(j,i));
          green=exp(%i.*k0.*(RffVx.*pos(1,:)+...
              RffVy.*pos(2,:)+RffVz.*pos(3,:)));
              
          OpxJx(1,:)= %i*k0/(4*%pi).*dS.*z0.*green.*(1-RffVx*RffVx);
          OpxJy(1,:)= %i*k0/(4*%pi).*dS.*z0.*green.*(-RffVy*RffVx);
          OpxJz(1,:)= %i*k0/(4*%pi).*dS.*z0.*green.*(-RffVz*RffVx);
          OpxMx(1,:)= zeros(1,size(pos,2));
          OpxMy(1,:)= %i*k0/(4*%pi).*dS.*green.*(RffVz);
          OpxMz(1,:)= %i*k0/(4*%pi).*dS.*green.*(-RffVy);
  
          OpyJx(1,:)= %i*k0/(4*%pi).*dS.*z0.*green.*(-RffVx*RffVy);
          OpyJy(1,:)= %i*k0/(4*%pi).*dS.*z0.*green.*(1-RffVy*RffVy);
          OpyJz(1,:)= %i*k0/(4*%pi).*dS.*z0.*green.*(-RffVz*RffVy);
          OpyMx(1,:)= %i*k0/(4*%pi).*dS.*green.*(-RffVz);
          OpyMy(1,:)= zeros(1,size(pos,2));
          OpyMz(1,:)= %i*k0/(4*%pi).*dS.*green.*(RffVx);
  
          OpzJx(1,:)= %i*k0/(4*%pi).*dS.*z0.*green.*(-RffVx*RffVz);
          OpzJy(1,:)= %i*k0/(4*%pi).*dS.*z0.*green.*(-RffVy*RffVz);
          OpzJz(1,:)= %i*k0/(4*%pi).*dS.*z0.*green.*(1-RffVz*RffVz);
          OpzMx(1,:)= %i*k0/(4*%pi).*dS.*green.*(RffVy);
          OpzMy(1,:)= %i*k0/(4*%pi).*dS.*green.*(-RffVx);
          OpzMz(1,:)= zeros(1,size(pos,2));
          
          Op.tJx(i,:,j)= cos(thetaFF(j,i))*cos(phiFF(j,i))*...
              OpxJx(1,:)+cos(thetaFF(j,i)).*sin(phiFF(j,i))*...
              OpyJx(1,:)-sin(thetaFF(j,i))*OpzJx(1,:);
          Op.tJy(i,:,j)= cos(thetaFF(j,i))*cos(phiFF(j,i))*...
              OpxJy(1,:)+cos(thetaFF(j,i)).*sin(phiFF(j,i))*...
              OpyJy(1,:)-sin(thetaFF(j,i))*OpzJy(1,:);
          Op.tJz(i,:,j)= cos(thetaFF(j,i))*cos(phiFF(j,i))*...
              OpxJz(1,:)+cos(thetaFF(j,i)).*sin(phiFF(j,i))*...
              OpyJz(1,:)-sin(thetaFF(j,i))*OpzJz(1,:);
          Op.tMx(i,:,j)= cos(thetaFF(j,i))*cos(phiFF(j,i))*...
              OpxMx(1,:)+cos(thetaFF(j,i)).*sin(phiFF(j,i))*...
              OpyMx(1,:)-sin(thetaFF(j,i))*OpzMx(1,:);
          Op.tMy(i,:,j)= cos(thetaFF(j,i))*cos(phiFF(j,i))*...
              OpxMy(1,:)+cos(thetaFF(j,i)).*sin(phiFF(j,i))*...
              OpyMy(1,:)-sin(thetaFF(j,i))*OpzMy(1,:);
          Op.tMz(i,:,j)= cos(thetaFF(j,i))*cos(phiFF(j,i))*...
              OpxMz(1,:)+cos(thetaFF(j,i)).*sin(phiFF(j,i))*...
              OpyMz(1,:)-sin(thetaFF(j,i))*OpzMz(1,:);
          
          Op.pJx(i,:,j)= - sin(phiFF(j,i))*OpxJx(1,:) + ...
              cos(phiFF(j,i))*OpyJx(1,:);
          Op.pJy(i,:,j)= - sin(phiFF(j,i))*OpxJy(1,:) + ...
              cos(phiFF(j,i))*OpyJy(1,:);
          Op.pJz(i,:,j)= - sin(phiFF(j,i))*OpxJz(1,:) + ...
              cos(phiFF(j,i))*OpyJz(1,:);
          Op.pMx(i,:,j)= - sin(phiFF(j,i))*OpxMx(1,:) + ...
              cos(phiFF(j,i))*OpyMx(1,:);
          Op.pMy(i,:,j)= - sin(phiFF(j,i))*OpxMy(1,:) + ...
              cos(phiFF(j,i))*OpyMy(1,:);
          Op.pMz(i,:,j)= - sin(phiFF(j,i))*OpxMz(1,:) + ...
              cos(phiFF(j,i))*OpyMz(1,:)
      end
  end
endfunction