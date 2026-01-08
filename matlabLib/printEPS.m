% Prints to PDF by EPSTOPDF call (better rendering than internal PDF 
% printer). Always print in vector format
%
% printEPS(path, filename, orientType)
%
% IN: path = path\to\directory
%     filename = string
% INopt: orientType = 'landscape' or 'portrait'
%
% Laurent Ntibarikure
function printEPS(path, filename, orientType)

figure(gcf);
if nargin>2
  orient(orientType);
end
set(gcf, 'Renderer', 'painters');
print('-depsc', [path, filename]);
system(['epstopdf ',filename, '.eps']);
system(['rm ',filename, '.eps']);