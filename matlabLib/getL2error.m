% Computes relative L2 (Euclidean) error norm
%
% DEFINITION:
%   Error_L2 = ||y_approx - y_ref||_2 / ||y_ref||_2
%            = sqrt(sum(|y_approx - y_ref|^2)) / sqrt(sum(|y_ref|^2))
%            [dimensionless ratio, typical: -20 to -100 dB]
%
% USE: Most common error metric. Amplitude-weighted, penalizes large errors.
%      Standard for radiation pattern comparison and SVD truncation studies.
%
% error = getL2error(matrix, refMatrix)
%
% IN: matrix = approximate/test data (vector or matrix)
%     refMatrix = reference/ground-truth data (same size)
%
% OUT: error = relative L2 error  [0 = perfect, 1 = completely wrong]
%              Printed in dB: 20*log10(error) for pattern metrics
%
% Laurent Ntibarikure
function error = getL2error(matrix, refMatrix)

error = norm(refMatrix-matrix)/norm(refMatrix);
fprintf('##> Relative L2 error = %2.4g\n', error);