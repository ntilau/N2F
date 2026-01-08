%%% DFT-truncation of the near to far field operator
clear; clc; close all;
addpath('../../matlabLib');

nbrElems_x = 9; % nbr of point sources for a linear array
N = 300; % nbr look angle samples - need oversampling to avoid aliasing !
Nbr = 11:4:31; % nbr of coefficients to retain in the DFT-truncation

arrayPos = buildArray(1, nbrElems_x, .5, 1, .5);
[xMin, xMax, yMin, yMax, zMin, zMax, xPts, yPts, zPts] = ...
  getBoxDim(1, arrayPos, .5, .1, .5 );
[boxPos, boxN, dS, mSize] = ...
  buildBox([1 1 1 1 1 1], xMin, xMax, yMin, yMax, zMin, zMax,...
  xPts, yPts, zPts, 1, 0, 0);
[Rmag, NdotRV] = getBoxVectors(arrayPos, boxPos, boxN);
excitPhasor = sf_excitations(1, arrayPos, 0, 0);
[psi, delPsi] = sf_nfSolver(1, excitPhasor, Rmag, NdotRV);
%% nf2ff operator and ff reference
thetaFF = linspace(0, 2*pi*(N-1)/N,N);
phiFF = 0;
[Lpsi, LdelPsi] = sf_nf2ffOperator(1, thetaFF, phiFF, ...
  boxPos, boxN, dS);
fPsiRef = sf_directffSolver(1, thetaFF, phiFF, ...
  excitPhasor, arrayPos);
%% DFT-truncation
Size = size(Lpsi);
LpsiTrunc = zeros(Size);
LdelPsiTrunc = zeros(Size);
Nbr = [Nbr,floor(Size(1)/2)+1];
Coeffs = zeros(1,length(Nbr));
TruncError = zeros(3,length(Nbr));
for n = 1:length(Nbr)
  disp('# FFT Truncation ...')
  Delta= Nbr(length(Nbr)+1-n);
  for i=1:size(Lpsi,3)
    %- operator on the near fields
    LpsiFFT = fft(Lpsi(:,:,i));
    LpsiFFT(Delta+2:Size(1)-Delta,:) = 0;
    LpsiTrunc(:,:,i) = ifft(LpsiFFT);
    clear LpsiFFT;
    %- operator on the near fields derivatives
    LdelPsiFFT = fft(LdelPsi(:,:,i));
    LdelPsiFFT(Delta+2:Size(1)-Delta,:) = 0;
    LdelPsiTrunc(:,:,i) = ifft(LdelPsiFFT);
    clear LdelPsiFFT;
  end

  %%-- Compute patterns
  approxPattern = zeros(size(fPsiRef));
  pattern = zeros(size(fPsiRef));
  for i=1:length(phiFF)
    approxPattern(i,:) =  LpsiTrunc(:,:,i) * psi.' + ...
      LdelPsiTrunc(:,:,i) * delPsi.';
    pattern(i,:) =  Lpsi(:,:,i) * psi.' +  LdelPsi(:,:,i) * delPsi.';
  end
  Coeffs(n) = 2*Delta+1;
  TruncError(1,n)= getL2error(approxPattern, pattern);
  TruncError(2,n)= getL2error(approxPattern, fPsiRef);
end
%% pattern L2Error
figure;
semilogy(Coeffs(2:length(Coeffs)), ...
  TruncError(1,2:length(TruncError(1,:))), '-*b','LineWidth', 1.5,...
      'MarkerSize', 7);
hold on;
semilogy(Coeffs(2:length(Coeffs)), ...
  TruncError(2,2:length(TruncError(2,:))), '-xk','LineWidth', 1.5,...
      'MarkerSize', 7);
xlabel('No of DFT coefficients', 'FontSize', 12)
ylabel('Relative Error', 'FontSize', 12)
legend('N2F', 'Direct');
axis('tight');
% printEPS('',['errorDFT',num2str(nbrElems_x)]);