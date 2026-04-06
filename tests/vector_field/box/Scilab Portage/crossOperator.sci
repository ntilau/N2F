function [Cx,Cy,Cz] = crossOperator(Ax,Ay,Az,Bx,By,Bz)
  Cx = Ay.*Bz - Az.*By;
  Cy = Az.*Bx - Ax.*Bz;
  Cz = Ax.*By - Ay.*Bx;
endfunction