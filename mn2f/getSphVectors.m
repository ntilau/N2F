% Computes the relevant vectors for near field to far field computations
% from a bounding sphere
%
% [N, rhonMag, NdotRhonV] = getSphVectors(arrayPos, spherePos)
%
% IN: arrayPos = Cartesian coordinates of the array elements
%     spherePos = Cartesian coordinates of the sphere sampling points
%
% OUT: rhonMag = vectors from each array element to the sphere sampling
%                points
%      NdotRhonV = for each array element, dot product between the normal
%                  unit vector of the surface patch and the unit vector
%                  pointing from the point source to the surface patch
%      N = normal unit vector to the patches, outwardly directed from the
%          sphere
%
% Laurent Ntibarikure
function [rhonMag, NdotRhonV, N] = getSphVectors(arrayPos, spherePos)
%% Normal versor of the patches in cartesian coordinates
radius = sqrt(spherePos(1,:).^2 + spherePos(2,:).^2 + spherePos(3,:).^2);
N(1,:) = spherePos(1,:) ./ radius;
N(2,:) = spherePos(2,:) ./ radius;
N(3,:) = spherePos(3,:) ./ radius;
%% Array related vectors to the patches
% For each elements vector to the integration surface, computation of the 
% distance to surface patches and the scalar product to their normal versor
rhonMag = zeros(size(arrayPos,2),size(spherePos,2));
NdotRhonV = zeros(size(arrayPos,2),size(spherePos,2));
for i=1:size(arrayPos,2)
  R(1,:) = spherePos(1,:) - arrayPos(1,i);
  R(2,:) = spherePos(2,:) - arrayPos(2,i);
  R(3,:) = spherePos(3,:);
  rhonMag(i,:) = sqrt( R(1,:).^2 + R(2,:).^2 + R(3,:).^2 );
  % scalar product of the versor to surface normal
  NdotRhonV(i,:) = dot(N,R)./ rhonMag(i,:);
end