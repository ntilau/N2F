//% Builds a box that encompass the array over which will be computed 
//% the near fields -> returns the nf sampling pts position and patch area
function [boxPos, boxN, dS] =...
    boxBuilder(k0, arrayPos, WLranging, WLspacing, planeMode, planeExt)
lambda0 = 2*%pi/k0;
planeExt = 2*planeExt*WLranging/WLspacing;
// Checking the choosen geometry
if(planeMode == 1)
    nbrFaces = 2;
else
    nbrFaces = 6;
end
// box geometry
if WLspacing > WLranging
    disp('Achtung!!! -->  ranging < spacing');
end
spacing = WLspacing*lambda0;
ranging = WLranging*lambda0;
arrayLength_x = max(arrayPos(1,:))-min(arrayPos(1,:));
arrayLength_y = max(arrayPos(2,:))-min(arrayPos(2,:));
dS = spacing^2;
// warning off all;
nbrSamples_x = floor((arrayLength_x+2*spacing)/spacing+1+(planeExt-1)*2 );
nbrSamples_y = floor((arrayLength_y+2*spacing)/spacing+1+(planeExt-1)*2 );
nbrSamples_z = floor(2*ranging/spacing+1);
// warning on all;
// Defining symmetries
if pmodulo(nbrSamples_x,2)
    halfNbrSamples_x = floor(nbrSamples_x/2);
else
    halfNbrSamples_x = (nbrSamples_x-1)/2;
end
if pmodulo(nbrSamples_y,2)
    halfNbrSamples_y = floor(nbrSamples_y/2);
else
    halfNbrSamples_y = (nbrSamples_y-1)/2;
end
if pmodulo(nbrSamples_z,2)
    halfNbrSamples_z = floor(nbrSamples_z/2);
else
    halfNbrSamples_z = (nbrSamples_z-1)/2;
end
// Building box
// XY face up
face(1).pos_x =  ones(nbrSamples_y,1) * ...
    (-halfNbrSamples_x:halfNbrSamples_x)* ...
    spacing;
face(1).pos_y = (-halfNbrSamples_y:halfNbrSamples_y)' * ...
    ones(1,nbrSamples_x) * ...
    spacing;
face(1).pos_z = (halfNbrSamples_z + 0.5) * ...
    spacing * ones(nbrSamples_y,nbrSamples_x);
face(1).nVx = 0;
face(1).nVy = 0;
face(1).nVz = 1;

// XY face down
face(2).pos_x =  face(1).pos_x;
face(2).pos_y =  face(1).pos_y;
face(2).pos_z =  - face(1).pos_z;
face(2).nVx = 0;
face(2).nVy = 0;
face(2).nVz = -1;

if(planeMode == 0)
    // XZ face up
    face(3).pos_x =  ones(nbrSamples_z,1) * ...
        (-halfNbrSamples_x:halfNbrSamples_x)* ...
        spacing;
    face(3).pos_y = (halfNbrSamples_y + 0.5) * ...
        spacing * ones(nbrSamples_z,nbrSamples_x);
    face(3).pos_z = ...
        (-halfNbrSamples_z:halfNbrSamples_z)' * ...
        ones(1,nbrSamples_x) * ...
        spacing;
    face(3).nVx = 0;
    face(3).nVy = 1;
    face(3).nVz = 0;

    // XZ face down
    face(4).pos_x =  face(3).pos_x;
    face(4).pos_y =  - face(3).pos_y;
    face(4).pos_z =  face(3).pos_z;
    face(4).nVx = 0;
    face(4).nVy = -1;
    face(4).nVz = 0;

    // YZ face up
    face(5).pos_x = (halfNbrSamples_x + 0.5) * ...
        spacing * ones(nbrSamples_z,nbrSamples_y);
    face(5).pos_y =  ones(nbrSamples_z,1) * ...
        (-halfNbrSamples_y:halfNbrSamples_y)* ...
        spacing;
    face(5).pos_z = (-halfNbrSamples_z:halfNbrSamples_z)' * ...
        ones(1,nbrSamples_y) * ...
        spacing;
    face(5).nVx = 1;
    face(5).nVy = 0;
    face(5).nVz = 0;

    // YZ face down
    face(6).pos_x =  - face(5).pos_x;
    face(6).pos_y =  face(5).pos_y;
    face(6).pos_z =  face(5).pos_z;
    face(6).nVx = -1;
    face(6).nVy = 0;
    face(6).nVz = 0;
end
//% Total number of samples
totalNbrSamples = 0;
for k=1:nbrFaces
    totalNbrSamples = totalNbrSamples + ...
        size(face(k).pos_x,1) * ...
        size(face(k).pos_x,2);
end
//%
boxPos = zeros(3,totalNbrSamples);
boxN = zeros(3,totalNbrSamples);
valTot = 0;
for k=1:nbrFaces
    val = size(face(k).pos_x(:),1);
    boxPos(1,valTot+(1:val))= (face(k).pos_x(:)).';
    boxPos(2,valTot+(1:val))= (face(k).pos_y(:)).';
    boxPos(3,valTot+(1:val))= (face(k).pos_z(:)).';
    boxN(1,valTot+(1:val))= face(k).nVx(:)*ones(1,val);
    boxN(2,valTot+(1:val))= face(k).nVy(:)*ones(1,val);
    boxN(3,valTot+(1:val))= face(k).nVz(:)*ones(1,val);
    valTot = valTot + val;
  end
  endfunction