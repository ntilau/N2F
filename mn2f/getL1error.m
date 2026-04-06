% Returns the relative error in the L1 norm between two vectors/matrices
%
% error = getL1error(matrix, refMatrix)
%
% IN: matrix = the vector/matrix for which we compute the error
%     refMatrix = reference matrix for relative error
%
% OUT: error = scalar value [non-dimensional]
%
% Laurent Ntibarikure
function error = getL1error(matrix, refMatrix)

error = norm(refMatrix-matrix,1)/norm(refMatrix,1);
fprintf('##> Relative L1 error = %2.4g\n', error);