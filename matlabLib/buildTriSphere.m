% Builds the triangulated sphere encompassing the array and the normal unit
% vectors to the triangles and outwardly directed from the array for near
% field to far field computation
%
% [triCent, tridS, triN] = buildTriSphere(lambda, ...
%   arrayPos, ext, triOrder, plotTriSphere)
%
% IN: radius = radius of the sphere in [m]
%     triOrder = triOrder of the triangulation
%     plotTriSphere = plots the sphere if asserted
%
% OUT: triCent = Cartesian components of the vectors pointing to the 
%                centroids of the triangles
%      tridS = areas of the triangles
%      triN = normals to the triangles
%
% Laurent Ntibarikure
function [triCent, tridS, triN] = buildTriSphere(radius, triOrder, ...
  plotTriSphere)
%% sphere structure construction
[p, t] = getTriSphMesh(triOrder);
p = (radius * eye(3) * p.').'; % scaling unit sphere radius to radius
nbrTri = size(t,1); % nbr of triangles
%% centroids & patches area & normal outwardly directed unit vectors
centx = zeros(1,nbrTri); centy = centx; centz = centx;
nVx = centx; nVy = centx; nVz = centx;
tridS = centx;
for i=1:nbrTri
  % each vertex is a column vector of x, y, z components    
  vert1 = (p(t(i,1),:)).';
  vert2 = (p(t(i,2),:)).';
  vert3 = (p(t(i,3),:)).';
  tri = [vert1,vert2,vert3];
  % centroids
  centx(i) = 1/3*sum(tri(1,:));
  centy(i) = 1/3*sum(tri(2,:));
  centz(i) = 1/3*sum(tri(3,:));
  % area of the triangle
  m1 = [tri(1,:);tri(2,:);ones(1,3)];
  m2 = [tri(2,:);tri(3,:);ones(1,3)];
  m3 = [tri(3,:);tri(1,:);ones(1,3)];
  tridS(i) = 1/2*sqrt( det(m1)^2 + det(m2)^2 + det(m3)^2 );
  % Normals
  v = vert3 - vert2;
  w = vert2 - vert1;
  nx = v(2)*w(3) - w(2)*v(3);
  ny = v(3)*w(1) - w(3)*v(1);
  nz = v(1)*w(2) - w(1)*v(2);
  nMag = sqrt(nx.^2 + ny.^2 + nz.^2);
  nVx(i) = nx ./ nMag;
  nVy(i) = ny ./ nMag;
  nVz(i) = nz ./ nMag;
end
triCent = [centx; centy; centz];
triN = [nVx; nVy; nVz];
%% plot sphere
if plotTriSphere 
  figure;
  trisurf(t, p(:,1), p(:,2), p(:,3), 'EdgeColor', [.1 .1 .3], ...
      'FaceColor', [.3 .3 .6] );
  axis equal; axis vis3d; view(3);
  xlabel('x');ylabel('y');zlabel('z');
end