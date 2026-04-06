% Exponential angle bisection strategy for adaptive N2F sampling
%
% Generates a sequence of angles optimized for low-rank SVD approximation of N2F operators.
% Uses exponential (power-of-2) bisection to concentrate sampling at critical angles.
%
% ALGORITHM: Binary tree angle refinement
%   Iteration 0: angles = [range/2, 0, range]  (3 points)
%   Iteration 1: insert midpoint between consecutive angles
%             Creates sequence like [0, c, range/2, c, range] (5 → 9 points)
%   Iteration 2: repeat bisection with new midpoints
%             Results in powers-of-2-spaced refinement: 1, 3, 9, 17, 33, 65, 129, 257...
%   Iteration N: 2^N + 1 angles (after sufficient iterations)
%
% RATIONALE:
%   - Exponential sequence [1 2 3 5 9 17 33 65 129 ...] arises naturally from bisection
%   - Pattern optimizes SVD low-rank models: concentrates DoF where field structure is richest
%   - Avoids regular linear sampling which may miss important spectral content
%   - Effective for pattern compression in radiation problems
%
% [angles, nbrAngles] = getSpanningAngles(maxNbrAngles, range)
%
% IN: maxNbrAngles = target truncation level (limits sequence length)
%     range = angular extent in degrees [0°, 180°], typically scan region around broadside
%            Example: range=90 → angles span primary lobe region
%
% OUT: angles = vector of spanning angle values in degrees (unsorted)
%      nbrAngles = cumulative count of angles at each bisection level
%                 Useful to track refinement progression
%
% EXAMPLE:
%   [angles, nbrAngles] = getSpanningAngles(40, 90)
%   angles =   [45 0 90 22.5 67.5 11.25 33.75 56.25 78.75 ...]
%   nbrAngles = [1  3   5    9   17   33   65   129]  (cum. at each iter)
%
% Laurent Ntibarikure
function [angles, nbrAngles] = getSpanningAngles(maxNbrAngles, range)
% Initialize with three strategic angles: range/2 (center), 0 (edge), range (edge)
angles = [range/2 0 range];
order = floor(log2(maxNbrAngles))-1;  % Determine bisection depth
dataInc = [1 2];  % Track cumulative angle count per iteration
nbrAngles = [1 3];  % Initial: 1 angle, then 3 after first iteration

for i=1:order
  % Sort current angles to enable bisection
  temp = sort(angles);
  % Find smallest spacing for next bisection
  new = temp(2)/2;  % Half of second-smallest angle creates new resolution
  angles(length(angles)+1) = new;  % Add new bisecting angle
  
  % Fill in all multiples of 'new' up to range
  data = 1;  % Count of newly added angles in this iteration
  Parts = range/new;  % Number of multiples to add
  for j=1:Parts
    % Only add if this angle not already present (avoid duplicates)
    if(~any(angles == new*j ))
      angles(length(angles)+1) = new*j;
      data = data + 1;
    end
  end
  dataInc(i+2) = data;  % Store iteration contributes 'data' new angles
  nbrAngles(i+2) = sum(dataInc(1:i+2));  % Cumulative total
end
