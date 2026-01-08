% plots the surface power density on the encompassing surface
%
% vf_plotSurfacePowerDensity(handle, mSize, surfPos, S)
%
% IN: handle = figure handle for plot update
%     mSize = matrix size for surface plot (surf())
%     surfPos = surface sampling points locations
%     S = Poynting vector on the surface
%
% Laurent Ntibarikure
function vf_plotSurfacePowerDensity(handle, mSize, surfPos, S)

figure(handle);
colormap(jet);

valTot = 0;
for k=1:size(mSize,3)
  val = mSize(1,1,k) * mSize(1,2,k);
  mS = vector2matrix(mSize(:,:,k), S(1,valTot+(1:val)).');
  pos_x = vector2matrix(mSize(:,:,k), surfPos(1,valTot+(1:val)).');
  pos_y = vector2matrix(mSize(:,:,k), surfPos(2,valTot+(1:val)).');
  pos_z = vector2matrix(mSize(:,:,k), surfPos(3,valTot+(1:val)).');
  surf(pos_x, pos_y, pos_z, 1/2*real(mS), 'FaceAlpha',.75, 'EdgeAlpha',.1, ...
    'HandleVisibility','off');
  hold on; axis('equal');
  valTot = valTot + val;
end