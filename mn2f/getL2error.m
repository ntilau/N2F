% Returns the relative error in the L2 norm between two vectors/matrices
%
% error = getL2error(matrix, refMatrix)
%
% IN: matrix = the vector/matrix for which we compute the error
%     refMatrix = reference matrix for relative error
%
% OUT: error = scalar value [non-dimensional]
%
% Laurent Ntibarikure
function error = getL2error(matrix, refMatrix)

error = norm(refMatrix-matrix)/norm(refMatrix);
fprintf('##> Relative L2 error = %2.4g\n', error);