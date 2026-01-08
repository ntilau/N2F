% Builds a box that encompass the array over which will be computed the
% near fields on the selected faces and returns the near field sampling
% points position and the patches areas
%
% [boxPos, boxN, dS, mSize] =...
%   buildBox(faces, xMin, xMax, yMin, yMax, zMin, zMax,...
%   xPts, yPts, zPts, scale, plotFlag, forPlot)
% 
% IN: faces = boolean vector of the faces to select [XY face up (top), 
%             XY face down (bottom), XZ face up, XZ face down,
%             YZ face up, YZ face down]
%     xMin, xMax, yMin, yMax, zMin, zMax = box dimensions
%     xPts, yPts, zPts = nbr of sampling points along x, y, z directions
%     scale = reduction of the box dimensions in [m] that scales
%             proportionally the box within the initial dimensions
%     plotFlag = plots the sampling grid if asserted
%     forPlot = allows proper box dimensions plot if asserted letting the
%               grid be related to the surface patches
%
% OUT: boxPos = positions of the sampling points or, grid intersections of
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