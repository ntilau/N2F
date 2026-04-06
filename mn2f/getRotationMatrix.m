% Computes the spherical rotation of "a" radiants of a vector x in
% cartesian coordinates to the vector x', considering the rotation axis
% defined by the cartesian components (c1,c2,c3)
%
% Rotation matrix =
%          |1 0 0|              |c1c1 c1c2 c1c3|        |  0 -c3  c2|
%     cos a|0 1 0| + (1 - cos a)|c2c1 c2c2 c2c3| + sin a| c3   0 -c1| 
%          |0 0 1|              |c3c1 c3c2 c3c3|        |-c2  c1   0|
%
% rotMatrix = getRotationMatrix(rotAngle,rotAxis)
%
% IN: rotAngle = rotation angle in degrees
%     rotAxis = [x y z] components of the unit vector pointing along the
%               rotation axis
%
% OUT: rotMatrix = 3x3 square matrix for linear transformation of vectors
%                  given in Cartesian components and restitutes rotated
%                  vectors
%
% Laurent Ntibarikure
function rotMatrix = getRotationMatrix(rotAngle,rotAxis)

rotAxis = rotAxis/norm(rotAxis);
rotAngle = deg2rad(rotAngle);
rotMatrix = cos(rotAngle) * eye(3) + ...
  (1-cos(rotAngle)) * rotAxis.' * ...
  rotAxis + sin(rotAngle) * ...
  [ 0, -rotAxis(3), rotAxis(2); ...
  rotAxis(3), 0 , -rotAxis(1); ...
  -rotAxis(2), rotAxis(1), 0 ];