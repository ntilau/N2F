//% array currents for beamforming
// the weighting parameter allows to apply windowing coefficients like
// Hamming and Kaiser ones to achieve a narrow beam
function [J,M] =  feedingCurrents(k0, z0, nbrElems_x, nbrElems_y, ...
   spacing_x, spacing_y, tJ, pJ, tM, pM, steering_x, steering_y, weighting)
exec('spherical2cartesian.sci',-1);
[i,j]=meshgrid(1:nbrElems_x,1:nbrElems_y);
steeringPhasor =  weighting .* ...  
    exp( - %i .* i .* k0 .* spacing_x .* sin(steering_x * %pi/180)) .* ...
    exp( - %i .* j .* k0 .* spacing_y .* sin(steering_y * %pi/180));
[Jx,Jy,Jz] = spherical2cartesian(steeringPhasor,0,0,tJ*%pi/180,pJ*%pi/180); 
[Mx,My,Mz] = spherical2cartesian(z0*steeringPhasor,0,0,tM*%pi/180,pM*%pi/180); 
J(:,1)= Jx(:); J(:,2)= Jy(:); J(:,3)= Jz(:); // [A/m]
M(:,1)= Mx(:); M(:,2)= My(:); M(:,3)= Mz(:); // [V/m]
endfunction