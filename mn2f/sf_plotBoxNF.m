% Plots the the near field distribuition (psi) on the box
%
% sf_plotBoxNF(mSize, boxPos, psi, chosenTitle)
%
% IN: mSize = to collect the near field vector into a matrix with
%             Ntheta rows and Nphi columns
%     boxPos = near field sampling points on the box
%     psi = vector of the near field sampled values
% INopt: chosenTitle = to put next to the |Psi| title
%
% Laurent Ntibarikure
function sf_plotBoxNF(mSize, boxPos, psi, chosenTitle)

nbrFaces = size(mSize,3);

valTot = 0;
for k=1:nbrFaces
  val = mSize(1,1,k) * mSize(1,2,k);
  mPsi = vector2matrix(mSize(:,:,k), psi(1,valTot+(1:val)).');
  pos_x = vector2matrix(mSize(:,:,k), boxPos(1,valTot+(1:val)).');
  pos_y = vector2matrix(mSize(:,:,k), boxPos(2,valTot+(1:val)).');
  pos_z = vector2matrix(mSize(:,:,k), boxPos(3,valTot+(1:val)).');
  surf(pos_x, pos_y, pos_z, abs(mPsi), 'FaceAlpha',.6,'EdgeAlpha',.3);
  hold on; axis('equal');
  valTot = valTot + val;
end

xlabel('x');
ylabel('y');
zlabel('z')
if nargin>3
  title(['|\Psi|', chosenTitle]);
end
axis tight;
view([30 45]);