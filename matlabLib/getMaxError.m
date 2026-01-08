% Returns the relative error in the L-inf norm between two vectors/matrices
%
% error = getMaxError(matrix, refMatrix)
%
% IN: matrix = the vector/matrix for which we compute the error
%     refMatrix = reference matrix for relative error
%
% OUT: error = scalar value [non-dimensional]
%
% Laurent Ntibarikure
function error = getMaxError(matrix, refMatrix)

error = norm(refMatrix-matrix, inf)/norm(refMatrix, inf);
fprintf('##> Relative max error = %2.4g\n', error);