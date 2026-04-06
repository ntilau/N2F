% Builds a triangulated sphere surface for near-field to far-field transformation
% using spherical harmonics triangulation with outward-pointing normals
%
% The triangulated sphere provides a discrete sampling surface for Huygens'
% principle or equivalent near-field computation with exact geometric properties:
% - Triangle centroids are sampling points (accurate field evaluation locations)
% - Triangle normals are outward-directed unit vectors needed for N2F integrals
% - Triangle areas enable proper numerical integration approximation
%
% GEOMETRIC COMPUTATION:
%   Triangle area: uses Heron's formula via determinants of projected components
%     Area = 0.5 * sqrt(||v1 x v2||²) = 0.5 * sqrt(det(xy)² + det(yz)² + det(zx)²)
%     This avoids direct cross product and is more numerically stable
%   Normal vector: n = (v2 - v1) x (v3 - v2) / ||...||  (right-hand rule for outward)
%
% [triCent, tridS, triN] = buildTriSphere(radius, triOrder, plotTriSphere)
%
% IN: radius = sphere radius [m]
%     triOrder = triangulation order/refinement level (affects mesh fineness)
%                Higher values: more triangles, finer mesh (e.g., 0,1,2,3,4)
%     plotTriSphere = boolean flag to visualize the generated triangulation
%
% OUT: triCent = Cartesian coordinates of triangle centroids (3 x nTriangles)
%                Centroid = (v1 + v2 + v3)/3 for each triangle
%      tridS = triangle areas [m²] (1 x nTriangles)
%              Used as integration weights in N2F operator computation
%      triN = outward unit normal vectors to triangles (3 x nTriangles)
%             Normalized: ||triN|| = 1.0, pointing away from sphere center
%
% DEPENDENCIES: getTriSphMesh(triOrder) - provides initial triangulation vertices/connectivity
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
  % Triangle area calculation using determinant projection method
  % Project 3D triangle onto XY, YZ, ZX planes and compute determinant of each 2D projection
  % m1 = [[x1 x2 x3]; [y1 y2 y3]; [1 1 1]] for XY-plane projection
  % Area = 0.5 * sqrt(det(m1)² + det(m2)² + det(m3)²)  [Heron's formula via Cayley-Menger]
  % This is more stable than direct cross product for small triangles
  m1 = [tri(1,:);tri(2,:);ones(1,3)];  % XY projection
  m2 = [tri(2,:);tri(3,:);ones(1,3)];  % YZ projection
  m3 = [tri(3,:);tri(1,:);ones(1,3)];  % ZX projection
  tridS(i) = 1/2*sqrt( det(m1)^2 + det(m2)^2 + det(m3)^2 );
  % Normal vector calculation using right-hand rule: n = (v3-v2) x (v2-v1)
  % Results in outward-pointing normal (toward positive sphere radius)
  % Magnitude must be normalized to 1.0 for proper far-field integration
  v = vert3 - vert2;  % edge vector
  w = vert2 - vert1;  % another edge vector
  % Cross product components: n = v x w
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