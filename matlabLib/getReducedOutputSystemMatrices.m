% Constructs reduced-order model (ROM) of N2F transformation operators
%
% OVERVIEW: Combines full electromagnetic system matrices (from FEM solver)
% with N2F operators to create compact output ROM for radiation patterns.
%
% WORKFLOW:
%   1. Load pre-computed FEM basis matrices from disk (comp2Ext, num2UnNum, etc.)
%   2. Assess coefficient importance: compute max||Opt[nbrCoeffs+1]|| to gauge
%      which DFT coefficients contribute significantly (nbrCoeffs levels checking)
%   3. Build rom matrices by chaining: Opt -> field functional -> FEM basis -> Q
%      Final ROM: rom = Opt * fieldFunctional * num2UnNum * comp2Ext * Q
%   4. Store as functions of observation angle phi (reduced dimensionality)
%
% INPUTS: Typical problem setup
%   System size from FEM: ~1000s unknowns → compressed to ~10s via ROM basis Q
%   Far-field samples: nTheta=256 angles × nPhi=4 phi cuts
%   DFT coefficients: nbrCoeffs=32 (typical, explores smooth patterns)
%
% OUTPUT: romCt, romCp can be applied as: E_ff(phi_k) = rom * field_snapshot
%         Enabling fast pattern sweeps without full FEM solve
%
% [romCt, romCp, romF, nbrSmpls, nbrCoeffs] = getReducedOutputSystemMatrices(...)
%
% IN: nbrCoeffs, nbrSmpls = DFT truncation parameters
%     k0, z0 = wavenumber and impedance
%     boxPos, boxN, dS = bounding box geometry parameters
%     phi = observation azimuth angles [rad]
%     Q = ROM basis matrix from prior reduced-order modeling
%
% OUT: romCt = reduced N2F operator for theta-polarized E-field (nCoeffs x nROM_basis x nPhi)
%      romCp = reduced N2F operator for phi-polarized E-field
%      romF, nbrSmpls, nbrCoeffs = diagnostic outputs
%
% DEPENDENCIES: Loads external FEM matrices from 'lte_fileset\' directory
%              Functions: n2fOpFieldsFFT(), mmread()
%
comp2Ext = mmread('lte_fileset\comp2Ext1.mm');
num2UnNum = mmread('lte_fileset\num2UnNum1.mm');
mFE = mmread('lte_fileset\functionalE.mm');
mFH = mmread('lte_fileset\functionalH.mm');
fieldFunctional  = [mFH;mFE];
% coeffLevels = 1;

fprintf('#> Computing the n2f Fields operators ... ')
tic();
% while minCoeffsLevel < coeffLevels
%   nbrSmpls = floor(nbrCoeffs * overSampling);
  coeffLevels = 0;
  for i=1:size(phi,2)
    [phiFF,thetaFF] = meshgrid(phi(1,i),...
      linspace(0,2*pi*(nbrSmpls-1)/nbrSmpls,nbrSmpls));
    [Opt, Opp] = n2fOpFieldsFFT(k0, z0, boxPos, boxN, dS, ...
        thetaFF, phiFF, nbrCoeffs);
    level = max(max(max(abs(Opt(nbrCoeffs+1,:,:)))), ...
      max(max(abs(Opp(nbrCoeffs+1,:,:)))));
    coeffLevels = max(coeffLevels,level);
  end
%   nbrCoeffs = nbrCoeffs + 1;
% end
fprintf('Elapsed %2.4g s.\n', toc());

% nbrCoeffs = nbrCoeffs - 1;
% nbrSmpls = floor(nbrCoeffs * overSampling);% floor(nbrCoeffs * overSampling);
fprintf('##> Coeff levels = %2.4g, nbrCoeffs = %d\n', ...
  coeffLevels, nbrCoeffs*2+1);

fprintf('#> Computing the Reduced Order n2f Fields operators ...\n')
tic();
romCt = zeros(nbrCoeffs*2+1,size(Q,2),size(phi,2));
romCp = romCt;
for i=1:size(phi,2)
  fprintf('##> Phi = %2.4g�\n', phi(1,i)*180/pi);
  [phiFF,thetaFF] = meshgrid(phi(1,i),...
    linspace(0,2*pi*(nbrSmpls-1)/nbrSmpls,nbrSmpls));
  [Opt, Opp] = n2fOpFieldsFFT(k0, z0, boxPos, boxN, dS, ...
      thetaFF, phiFF, nbrCoeffs);
  romCt(:,:,i) = Opt*fieldFunctional*num2UnNum*comp2Ext*Q;
  romCp(:,:,i) = Opp*fieldFunctional*num2UnNum*comp2Ext*Q;
  
end
romF = fieldFunctional*num2UnNum*comp2Ext*Q;
fprintf('Elapsed %2.4g s.\n', toc());
