% theta exponential angles' selection
% builds the sequence 1 2 3 5 9 17 33 65 129 ...
%
% [angles, nbrAngles] = getSpanningAngles(maxNbrAngles, range)
%
% IN: maxNbrAngles = truncation of the sequence
%     range = extention around broadside in degrees in [0°,180°]
%
% OUT: angles = theta values in degrees
%      nbrAngles = number of angles effectively computed
%
% Laurent Ntibarikure
function [angles, nbrAngles] = getSpanningAngles(maxNbrAngles, range)
angles = [range/2 0 range];
order = floor(log2(maxNbrAngles))-1;
dataInc = [1 2];
nbrAngles = [1 3];
for i=1:order
  temp = sort(angles);
  new = temp(2)/2;
  angles(length(angles)+1) = new;
  data = 1;
  Parts = range/new;
  for j=1:Parts
    if(~any(angles == new*j )) 
      angles(length(angles)+1) = new*j;
      data = data + 1; 
    end
  end
  dataInc(i+2) = data ;
  nbrAngles(i+2) = sum(dataInc(1:i+2));
end
