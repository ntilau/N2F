% Prints to PDF with internal PDF printer. Let Matlab choose raster or
% vector format for printing
%
% printPDF(path, filename, orientType)
%
% IN: path = path\to\directory
%     filename = string
% INopt: orientType = 'landscape' or 'portrait'
%
% Laurent Ntibarikure
function printPDF(path, filename, orientType)

figure(gcf);
if nargin>2
  orient(orientType);
end
print('-dpdf', [path, filename]);