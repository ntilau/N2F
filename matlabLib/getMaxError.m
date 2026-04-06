% Computes relative maximum (Chebyshev/L-infinity) error norm
%
% DEFINITION:
%   Error_max = ||y_approx - y_ref||_inf / ||y_ref||_inf
%             = max(|y_approx - y_ref|) / max(|y_ref|)  [dimensionless ratio]
%
% USE: Worst-case error analysis. Most stringent metric - single outlier dominates.
%      Critical for sidelobe control, null steering, and minimax design.
%
% error = getMaxError(matrix, refMatrix)
%
% IN: matrix = approximate/test data (vector or matrix)
%     refMatrix = reference/ground-truth data (same size)
%
% OUT: error = relative maximum error  [indicates worst-point deviation]
%             Single largest per-element error normalized
%
% Laurent Ntibarikure
function error = getMaxError(matrix, refMatrix)

error = norm(refMatrix-matrix, inf)/norm(refMatrix, inf);
fprintf('##> Relative max error = %2.4g\n', error);