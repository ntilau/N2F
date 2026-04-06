% Error induced to the pattern and singular values plot after low rank 
% approximation of the near field
%
% plotSVDerror(sPsi, nbrVectors, nError, fError, fRefError, ...
%   nbrElems_x, plotErrorOnly)
%
% IN: sPsi = diagonal matrix of the singular values
%     nbrVectors = number of near field pictures collected
%     nError = L2 error in the near field relatively to non tested
%     fError = L2 error in the transformed far field relatively to non
%              tested near fields
%     fRefError = L2 error in the transformed far field relatively to
%                 directly computed far field
%     nbrElems_x = number of point sources in the x direction (scanning
%                  direction)
%     plotErrorOnly = avoid plotting the singular values
%
% Laurent Ntibarikure 
 function plotSVDerror(sPsi, nbrVectors, nError, fError, fRefError, ...
  nbrElems_x, plotErrorOnly)

plot1 = false;
if nargin > 5
  if plotErrorOnly
    plot1 = true;
  end
end

nError = nError(nbrVectors);
fError = fError(nbrVectors);

figProp = getFigureProperties();
figure();
if ~plot1
  subplot(1,2,1);
end
semilogy(nbrVectors, nError,'+-b', 'LineWidth', figProp.lw, ...
  'MarkerSize', figProp.ms);
hold on
semilogy(nbrVectors, fError,'*-r', 'LineWidth', figProp.lw, ...
  'MarkerSize', figProp.ms);
maxError = max(nError(length(nbrVectors)),fError(length(nbrVectors)));
minError = min(nError(length(nbrVectors)),fError(length(nbrVectors)));
if nargin > 4
  fRefError = fRefError(nbrVectors);
  semilogy(nbrVectors, fRefError,'x-k', 'LineWidth', figProp.lw, ...
    'MarkerSize', figProp.ms);
end
axis('tight');
v = axis;
axis([v(1) v(2) 1e-15 1e5]);
if 10^(log10(minError)+(log10(maxError/minError)/2)) > ...
    10^(log10(v(3))+(log10(v(4)/v(3))/2))
  location = 'SouthEast';
else
  location = 'NorthEast';
end
xlabel('No of vectors', 'FontSize', figProp.fs);
ylabel('Relative error', 'FontSize', figProp.fs);
if nargin > 4
  legend('Near field', 'N2F far field', 'Direct far field', ...
    'Location',location);
else
  legend('Near field', 'N2F far field','Location',location);
end
if ~plot1
  subplot(1,2,2);
  semilogy(diag(sPsi),'o-b', 'LineWidth', figProp.lw, ...
    'MarkerSize', figProp.ms);
  hold on;
  if nargin > 5
    if nbrVectors(length(nbrVectors))> nbrElems_x
      text(nbrElems_x+1,sPsi(nbrElems_x,nbrElems_x),[...
      '\sigma_{',num2str(nbrElems_x), ...
      '} = ', sprintf('%2.4g',sPsi(nbrElems_x,nbrElems_x))]);
    else
      text(2,sPsi(nbrElems_x,nbrElems_x),[...
      '\sigma_{',num2str(nbrElems_x), ...
      '} = ', sprintf('%2.4g',sPsi(nbrElems_x,nbrElems_x))]);
    end
  end
  axis tight;
  v = axis;
  axis([v(1) v(2) 1e-15 1e5]);
  xlabel('n', 'FontSize', figProp.fs);
  ylabel('Singular values \sigma_n', 'FontSize', figProp.fs);
end