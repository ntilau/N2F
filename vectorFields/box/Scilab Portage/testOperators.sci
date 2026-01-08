// test of the operators validity
clear; clc;
//gstacksize('max');
exec('arrayBuilder.sci',-1);
exec('boxBuilder.sci',-1);
exec('nearFieldSolver.sci',-1);
exec('nf2equivCurrents.sci',-1);
exec('farFieldSolver.sci',-1);
exec('plotFF.sci',-1);
tic();
// electrical parameters
eps0=8.85418781761e-12;
mu0=%pi*4e-7;
c0=1/sqrt(eps0*mu0);
frequency=1e9;
z0=sqrt(mu0/eps0); // free-space impedance
k0=2*%pi*frequency/c0; // wavenumber
clear eps0 mu0 c0 frequency;
// array
weighting = 1; // fixed feeding current phasors [nbrElems_x,nbrElems_y]
[arrayPos, J, M] = arrayBuilder(k0, z0, 11, .5, 7, .5, 90, 270, 90, 0, ...
  30, 0, weighting); // Huygens' sources
// box
[boxPos, boxN, dS] =  boxBuilder(k0, arrayPos, .5, .1, 0, 2);
// nf
[boxEt, boxHt] = nearFieldSolver(k0, z0, arrayPos, boxPos, J, M);
clear arrayPos J M ;
[Js, Ms, Pr] = nf2equivCurrents(boxEt, boxHt, boxN, dS);
clear boxEt boxHt boxN;
printf('##> nearfields computation performed in %g s.\n', toc());
tic();
[gaint, gainp, thetaFF, phiFF] = farFieldSolver(k0, z0, boxPos, dS, .5, Js, Ms, Pr);
printf('##> farfields computation performed in %g s.\n', toc());


plotFF(thetaFF, 80, 1, gaint, gainp);