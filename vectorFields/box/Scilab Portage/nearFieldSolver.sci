//% Fields computation over the box
function [boxEt, boxHt] = nearFieldSolver(k0, z0, arrayPos, boxPos, J, M)
  
  exec('crossOperator.sci',-1);
  boxEt=zeros(3,size(boxPos,2));
  boxHt=zeros(3,size(boxPos,2));
  for i=1:size(arrayPos,2)
      Rx = boxPos(1,:) - arrayPos(1,i);
      Ry = boxPos(2,:) - arrayPos(2,i);
      Rz = boxPos(3,:);
      R=sqrt(Rx.^2+Ry.^2+Rz.^2);
      RxV=Rx./R; RyV=Ry./R; RzV=Rz./R;
      JdotRV=J(i,1) .* RxV+J(i,2) .* RyV+J(i,3) .* RzV;
      [JcrossRVx,JcrossRVy,JcrossRVz]= ...
          crossOperator(J(i,1),J(i,2),J(i,3),RxV,RyV,RzV);
      MdotRV=M(i,1) .* RxV+M(i,2) .* RyV+M(i,3) .* RzV;
      [McrossRVx,McrossRVy,McrossRVz]= ... 
      crossOperator(M(i,1),M(i,2),M(i,3),RxV,RyV,RzV);
      
      
      boxEt(1,:) = boxEt(1,:) + (z0*k0^2/(4*%pi) .* ( J(i,1) .* ...
          (-%i./(k0*R) - (k0*R).^-2 + %i.*(k0*R).^-3) + ...
          JdotRV .* RxV .* (%i./(k0*R)+3.*(k0*R).^-2-3*%i.*(k0*R).^-3) ) - ...
          k0^2/(4*%pi) .*  McrossRVx .* (%i./(k0 .* R)+1.*(k0 .* R).^-2) ) .* ...
          exp(-%i*k0 .* R);
      boxEt(2,:) = boxEt(2,:) + (z0*k0^2/(4*%pi) .* ( J(i,2) .* ...
          (-%i./(k0*R) - (k0*R).^-2 + %i.*(k0*R).^-3) + ...
          JdotRV .* RyV .* (%i./(k0*R)+3.*(k0*R).^-2-3*%i.*(k0*R).^-3) ) - ...
          k0^2/(4*%pi) .* McrossRVy .* (%i./(k0 .* R)+1.*(k0 .* R).^-2) ) .* ...
          exp(-%i*k0 .* R);
      boxEt(3,:) = boxEt(3,:) + ...
          ( z0*k0^2/(4*%pi) .* ( J(i,3) .* (-%i./(k0*R) - (k0*R).^-2 + %i.*(k0*R).^-3) +...
          JdotRV .* RzV .* (%i./(k0*R)+3.*(k0*R).^-2-3*%i.*(k0*R).^-3) )-...
          k0^2/(4*%pi) .* McrossRVz .* (%i./(k0 .* R)+1.*(k0 .* R).^-2) )  .* ...
          exp(-%i*k0 .* R);
      boxHt(1,:) = boxHt(1,:) + ...
          ( k0^2/(4*%pi*z0) .* ( M(i,1) .* (-%i./(k0*R) - (k0*R).^-2 + %i.*(k0*R).^-3) +...
          MdotRV .* RxV .* (%i./(k0*R)+3.*(k0*R).^-2-3*%i.*(k0*R).^-3) ) + ...
          k0^2/(4*%pi) .* JcrossRVx .* (%i./(k0 .* R)+1.*(k0 .* R).^-2) )   .* ...
          exp(-%i*k0 .* R);
      boxHt(2,:) = boxHt(2,:) + ...
          ( k0^2/(4*%pi*z0) .* ( M(i,2) .* (-%i./(k0*R) - (k0*R).^-2 + %i.*(k0*R).^-3) +...
          MdotRV .* RyV .* (%i./(k0*R)+3.*(k0*R).^-2-3*%i.*(k0*R).^-3) ) + ...
          k0^2/(4*%pi) .*  JcrossRVy .* (%i./(k0 .* R)+1.*(k0 .* R).^-2) ) .* ...
          exp(-%i*k0 .* R);
      boxHt(3,:) = boxHt(3,:) + ...
          (k0^2/(4*%pi*z0) .* ( M(i,3) .* (-%i./(k0*R) - (k0*R).^-2 + %i.*(k0*R).^-3) + ...
          MdotRV .* RzV .* (%i./(k0*R)+3.*(k0*R).^-2-3*%i.*(k0*R).^-3) ) + ...
          k0^2/(4*%pi) .*  JcrossRVz .* (%i./(k0 .* R)+1.*(k0 .* R).^-2) )  .* ...
          exp(-%i*k0 .* R);
  end
  
endfunction
  