% Converts a column vector back to original matrix form via linear indexing
%
% Inverse operation of MATLAB's matrix(:) syntax, which flattens by columns.
% Used to reconstruct 2D operator matrices from stacked N2F computation results.
%
% LINEAR INDEXING: MATLAB stores matrices in column-major order
%   For dim = [3, 4]: linear_index = row + (col-1)*nRows
%   Example: matrix(2,3) has linear_index = 2 + (3-1)*3 = 8
%
% CONVERSION: Given flattened vector form (from matrix(:)):
%   Extract column j: matrix(:,j) = vector(1:dim(1) + (j-1)*dim(1))
%   This retrieves contiguous block for each column
%
% matrix =  vector2matrix(dim, vector)
%
% IN: dim = [nRows, nCols] - target matrix dimensions
%     vector = column vector from matrix(:) flattening (length = nRows*nCols)
%
% OUT: matrix = reconstructed matrix of size (nRows x nCols)
%               Satisfies: vector = matrix(:) [identity after conversion]
%
% NOTE: Could use matrix = reshape(vector, dim) for cleaner MATLAB code
%       Current implementation preserves original architecture compatibility
%
% Laurent Ntibarikure
function matrix =  vector2matrix(dim, vector)
% Reconstruct matrix from column-vector via linear indexing of MATLAB's column-major storage
matrix = zeros(dim);
for j=1:dim(2)  % Iterate over columns
  % Extract rows for column j: indices = [start:start+nRows-1]
  % Linear indexing formula: 1-indexed row + (col-1)*nRows
  start_idx = (j-1)*dim(1) + 1;  % Starting linear index for column j
  end_idx = j*dim(1);            % Ending linear index for column j
  matrix(:,j) = vector(start_idx:end_idx,1);
end