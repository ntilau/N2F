//% Builds the planar array on the XY plane and computes the feeding
//% currents
function [arrayPos,J,M] = arrayBuilder(k0, z0, nbrElems_x, WLspacing_x, nbrElems_y, WLspacing_y, tJ, pJ, tM, pM, steering_x, steering_y, weighting)
exec('feedingCurrents.sci',-1);
lambda0 = 2*%pi/k0;
// array main geometry
spacing_x = WLspacing_x*lambda0;
spacing_y = WLspacing_y*lambda0;
// Symmetry checking
if pmodulo(nbrElems_x,2)
    halfNbrElems_x = floor(nbrElems_x/2);
else
    halfNbrElems_x = (nbrElems_x-1)/2;
end
if pmodulo(nbrElems_y,2)
    halfNbrElems_y = floor(nbrElems_y/2);
else
    halfNbrElems_y = (nbrElems_y-1)/2;
end 
// Determine array elements' position
pos_x = (-halfNbrElems_x:halfNbrElems_x)*spacing_x;
pos_y = (-halfNbrElems_y:halfNbrElems_y)*spacing_y;
[pos_x,pos_y] = meshgrid(pos_x,pos_y);
arrayPos(1,:)= pos_x(:).';
arrayPos(2,:)= pos_y(:).';
// array feeding currents
[J,M] =  feedingCurrents(k0, z0, nbrElems_x, nbrElems_y, spacing_x, spacing_y, tJ, pJ, tM, pM, steering_x, steering_y, weighting);
endfunction