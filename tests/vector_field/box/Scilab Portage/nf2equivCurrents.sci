function [Js, Ms, Pr] = nf2equivCurrents(boxEt, boxHt, boxN, dS)

Js=zeros(boxN);
Ms=zeros(boxN);
S=zeros(boxN);
exec('crossOperator.sci',-1);

  [x,y,z] = crossOperator(boxN(1,:),boxN(2,:),boxN(3,:), ...
      boxHt(1,:),boxHt(2,:),boxHt(3,:));
  Js(1,:)=x;
  Js(2,:)=y;
  Js(3,:)=z;
  [x,y,z] = crossOperator(boxEt(1,:),boxEt(2,:),boxEt(3,:),...
      boxN(1,:),boxN(2,:),boxN(3,:));
  Ms(1,:)=x;
  Ms(2,:)=y;
  Ms(3,:)=z;
  [x,y,z] = crossOperator(Js(1,:),Js(2,:),Js(3,:),...
      conj(Ms(1,:)),conj(Ms(2,:)),conj(Ms(3,:)));
  S(1,:)=x;
  S(2,:)=y;
  S(3,:)=z;
  Sr = sum(S.*boxN,'r');
  Pr=1/2*real(sum(Sr))*dS;
endfunction