% Constructs a rectangular bounding box surface for near-field sampling
%
% Creates rectangular faces with optional selective inclusion for N2F transformations.
% Enables flexible bounding surface choices (full box, partial enclosure, single plane)
%
% FACE INDEXING: faces = [faceTop, faceBottom, faceUp, faceDown, faceLeft, faceRight]
%   All indices use 1 = include, 0 = exclude
%   faceTop = XY plane at z=zMax (outward normal +z)
%   faceBottom = XY plane at z=zMin (outward normal -z)
%   faceUp = XZ plane at y=yMax (+y normal)
%   faceDown = XZ plane at y=yMin (-y normal)
%   faceLeft = YZ plane at x=xMin (-x normal)
%   faceRight = YZ plane at x=xMax (+x normal)
%
%   Example: [1 1 1 1 1 1] = all six faces (closed box)
%   Example: [0 0 1 1 1 1] = four sides only (open top/bottom)
%   Example: [1 0 0 0 0 0] = single face (for MoM problems)
%
% [boxPos, boxN, dS, mSize] = buildBox(faces, xMin, xMax, yMin, yMax, zMin, zMax, ...
%   xPts, yPts, zPts, scale, plotFlag, forPlot)
%
% IN: faces = [6x1] binary flag vector selecting faces as described above
%     xMin, xMax, yMin, yMax, zMin, zMax = box boundaries [m]
%     xPts, yPts, zPts = number of sampling points along each edge direction
%     scale = box contraction factor [m]; final dimensions = original - scale
%             Used for positioning inside room or near geometry
%     plotFlag = 1 to visualize the sampling grid
%     forPlot = 1 to ensure patches align for proper 3D visualization
%
% OUT: boxPos = Cartesian coordinates of sampling points (3 x nTotalPts)
%      boxN = outward normal unit vectors to each point (3 x nTotalPts)
%      dS = area per patch element [m^2] (1 x nTotalPts)
%           Uniform spacing: dS(i) = (dX * dY) or (dX * dZ) or (dY * dZ)
%      mSize = mesh metadata: dimensions, statistics, face tallies
%
% GEOMETRY: Each face is a regular rectangular mesh
%   XY face: dX = (xMax-xMin)/(xPts-1), dY = (yMax-yMin)/(yPts-1)
%   Outward normals always point AWAY from box interior
%
% Laurent Ntibarikure
%               the sampling patches if forPlot is asserted
%      boxN = outwardly directed normal unit vectors to the faces
%      dS = surface patches areas
%      mSize = matrices sizes collected in a 3D matrix of depth the number
%              of faces selected
%
% Laurent Ntibarikure 
function [boxPos, boxN, dS, mSize] =...
    buildBox(faces, xMin, xMax, yMin, yMax, zMin, zMax,...
    xPts, yPts, zPts, scale, plotFlag, forPlot)

nbrFaces = length(find(faces));
 
xVect = linspace(xMin, xMax, xPts+1);
yVect = linspace(yMin, yMax, yPts+1);
zVect = linspace(zMin, zMax, zPts+1);

dx = abs(xVect(1)-xVect(2));
dy = abs(yVect(1)-yVect(2));
dz = abs(zVect(1)-zVect(2));

if nargin>12 && forPlot
else
  xMin = xMin * scale;
  yMin = yMin * scale;
  zMin = zMin * scale;
  xMax = xMax * scale;
  yMax = yMax * scale;
  zMax = zMax * scale;
  
  xMinTmp = xMin + dx/2;
  yMinTmp = yMin + dy/2;
  zMinTmp = zMin + dz/2;
  xMaxTmp = xMax - dx/2;
  yMaxTmp = yMax - dy/2;
  zMaxTmp = zMax - dz/2;

  xVect = linspace(xMinTmp, xMaxTmp, xPts);
  yVect = linspace(yMinTmp, yMaxTmp, yPts);
  zVect = linspace(zMinTmp, zMaxTmp, zPts);
end
%% Building box
k=0;
if faces(1)
  % XY face up
  k=k+1;
  [face(k).pos_x, face(k).pos_y] = meshgrid(xVect, yVect);
  face(k).pos_z = (zMax)*ones(size(face(k).pos_x));
  face(k).dS = dx*dy*ones(size(face(k).pos_x));
  face(k).nVx = 0;
  face(k).nVy = 0;
  face(k).nVz = 1;
end

if faces(2)
  % XY face down
  k=k+1;
  [face(k).pos_x, face(k).pos_y] = meshgrid(xVect, yVect);
  face(k).pos_z = (zMin)*ones(size(face(k).pos_x));
  face(k).dS = dx*dy*ones(size(face(k).pos_x));
  face(k).nVx = 0;
  face(k).nVy = 0;
  face(k).nVz = -1;
end

if faces(3)
  % XZ face up
  k=k+1;
  [face(k).pos_x, face(k).pos_z] = meshgrid(xVect, zVect);
  face(k).pos_y = (yMax)*ones(size(face(k).pos_x));
  face(k).dS = dx*dz*ones(size(face(k).pos_x));
  face(k).nVx = 0;
  face(k).nVy = 1;
  face(k).nVz = 0;
end

if faces(4)
  % XZ face down
  k=k+1;
  [face(k).pos_x, face(k).pos_z] = meshgrid(xVect, zVect);
  face(k).pos_y = (yMin)*ones(size(face(k).pos_x));
  face(k).dS = dx*dz*ones(size(face(k).pos_x));
  face(k).nVx = 0;
  face(k).nVy = -1;
  face(k).nVz = 0;
end

if faces(5)
  % YZ face up
  k=k+1;
  [face(k).pos_y, face(k).pos_z] = meshgrid(yVect, zVect);
  face(k).pos_x = (xMax)*ones(size(face(k).pos_y));
  face(k).dS = dy*dz*ones(size(face(k).pos_y));
  face(k).nVx = 1;
  face(k).nVy = 0;
  face(k).nVz = 0;
end

if faces(6)
  % YZ face down
  k=k+1;
  [face(k).pos_y, face(k).pos_z] = meshgrid(yVect, zVect);
  face(k).pos_x = (xMin)*ones(size(face(k).pos_y));
  face(k).dS = dy*dz*ones(size(face(k).pos_y));
  face(k).nVx = -1;
  face(k).nVy = 0;
  face(k).nVz = 0;
end
%% Total number of samples
totalNbrSamples = 0;
if plotFlag
  figure;
  colormap(getColorMap);
end
for k=1:nbrFaces
  totalNbrSamples = totalNbrSamples + ...
    size(face(k).pos_x,1) * ...
    size(face(k).pos_x,2);
  if plotFlag
    mesh(face(k).pos_x,face(k).pos_y,face(k).pos_z, ...
      'FaceAlpha',0,'EdgeAlpha',1, 'EdgeColor', getColorMap);
    hold on; axis('equal');
  end
end
%%
boxPos = zeros(3,totalNbrSamples);
dS = zeros(1,totalNbrSamples);
boxN = zeros(3,totalNbrSamples);
mSize = zeros(1,2,nbrFaces);
valTot = 0;
for k=1:nbrFaces
  val = size(face(k).pos_x(:),1);
  mSize(:,:,k) = size(face(k).pos_x);
  boxPos(1,valTot+(1:val))= (face(k).pos_x(:)).';
  boxPos(2,valTot+(1:val))= (face(k).pos_y(:)).';
  boxPos(3,valTot+(1:val))= (face(k).pos_z(:)).';
  dS(1,valTot+(1:val))= (face(k).dS(:)).';
  boxN(1,valTot+(1:val))= face(k).nVx(:)*ones(1,val);
  boxN(2,valTot+(1:val))= face(k).nVy(:)*ones(1,val);
  boxN(3,valTot+(1:val))= face(k).nVz(:)*ones(1,val);
  valTot = valTot + val;
end