% Computes relative L1 (Manhattan/taxicab) error norm
%
% DEFINITION:
%   Error_L1 = ||y_approx - y_ref||_1 / ||y_ref||_1
%            = sum(|y_approx - y_ref|) / sum(|y_ref|)  [dimensionless ratio]
%
% USE: Amplitude/magnitude errors, robust to outliers, used for pattern fidelity
%      Less sensitive to phase errors than L2 norm
%
% error = getL1error(matrix, refMatrix)
%
% IN: matrix = approximate/test data (vector or matrix)
%     refMatrix = reference/ground-truth data (same size)
%
% OUT: error = relative L1 error  [0 = perfect, 1 = completely wrong]
%
% Laurent Ntibarikure
function error = getL1error(matrix, refMatrix)

error = norm(refMatrix-matrix,1)/norm(refMatrix,1);
fprintf('##> Relative L1 error = %2.4g\n', error);