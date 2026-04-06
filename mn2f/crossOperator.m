% Cross product of vectorial fields A and B
%
% [Cx,Cy,Cz] = crossOperator(Ax,Ay,Az,Bx,By,Bz)
%
% IN: Ax, Ay, Az = scalar or matrix of components (vector field)
%     Bx, By, Bz = scalar or matrix of components (vector field)
%
% OUT: Cx, Cy, Cz = dot product components of C = A . B
%
% Note: - can combine scalars to matrices: e.g. A is a single vector and B
%       is a set of vectors
%       - Do not care of the size of the matrices as the components are
%       already split
%
% Laurent Ntibarikure
function [Cx,Cy,Cz] = crossOperator(Ax,Ay,Az,Bx,By,Bz)
Cx = Ay.*Bz - Az.*By;
Cy = Az.*Bx - Ax.*Bz;
Cz = Ax.*By - Ay.*Bx;