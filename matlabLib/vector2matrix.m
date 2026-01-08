% Vector conversion to matrix. The vector is previouly obtained by the
% matrix we need to reconstruct with the command vector = matrix(:)
%
% matrix =  vector2matrix(dim, vector)
%
% IN: dim = (rows, columns)
%     vector = input data to transform into matrix
%
% OUT: matrix = vector data collected into a matrix of size "dim"
%
% Laurent Ntibarikure
function matrix =  vector2matrix(dim, vector)

% fprintf('-> Vector 2 Matrix conversion...\t')
% tic();
matrix = zeros(dim);
for j=1:dim(2)                
    matrix(:,j) = vector((1:dim(1))+(j-1)*dim(1),1);
end
% fprintf('Elapsed %2.4g s\n', toc());