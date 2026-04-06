% Computes the relevant vectors for scalar near field computations
%
% [rhonMag, boxNdotRhonV] = getBoxVectors(arrayPos, boxPos, boxN)
%
% IN: arrayPos = array elements positions
%     boxPos = box sampling points positions
%     boxN = outwardly directed normal unit vectors to the faces
%
% OUT: rhonMag = vector magnitude from the array positions to the box near
%                field sampling points
%      boxNdotRhonV = dot product between the normal unit vectors and the
%                     unit vectors in the near field sampling directions
%
% Laurent Ntibarikure
function [rhonMag, boxNdotRhonV] = getBoxVectors(arrayPos, boxPos, boxN)
%% Array related vectors to the patches
% For each elements vector to the integration surface, computation of the 
% distance to surface patches and the scalar product to their normal versor
rhonMag = zeros(size(arrayPos,2),size(boxPos,2));
boxNdotRhonV = zeros(size(arrayPos,2),size(boxPos,2));
for i=1:size(arrayPos,2)
  R(1,:) = boxPos(1,:) - arrayPos(1,i);
  R(2,:) = boxPos(2,:) - arrayPos(2,i);
  R(3,:) = boxPos(3,:);
  rhonMag(i,:) = sqrt( R(1,:).^2 + R(2,:).^2 + R(3,:).^2 );
  % scalar product of the versor to surface normal
  boxNdotRhonV(i,:) = dot(boxN,R)./ rhonMag(i,:);
end