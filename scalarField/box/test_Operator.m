%%% low-rank approximation of the N2F operator for a bounding sphere 
clear; clc; close all;
addpath('../../matlabLib');

arrayPos = buildArray(1, 9, .5, 1, .5);
[xMin, xMax, yMin, yMax, zMin, zMax, xPts, yPts, zPts] = ...
  getBoxDim(1, arrayPos, .5, .1, .5 );
[boxPos, boxN, dS, mSize] = ...
  buildBox([1 1 1 1 1 1], xMin, xMax, yMin, yMax, zMin, zMax,...
  xPts, yPts, zPts, 1, 0, 0);
[Rmag, NdotRV] = getBoxVectors(arrayPos, boxPos, boxN);
excitPhasor = sf_excitations(1, arrayPos, 0, 0);
[psi, delPsi] = sf_nfSolver(1, excitPhasor, Rmag, NdotRV);
%% nf2ff
dthetaFF = 1; % ff pattern resolution [°]
N = 300;
thetaFF = linspace(0, 2*pi*(N-1)/N,N);
phiFF = deg2rad([0 90]);
[Lpsi, LdelPsi] = sf_nf2ffOperator(1, thetaFF, phiFF, ...
  boxPos, boxN, dS);
fPsiRef = sf_directffSolver(1, thetaFF, phiFF, ...
  excitPhasor, arrayPos);
%% ----- plots
fPsi = zeros(length(phiFF),length(thetaFF));
for i=1:length(phiFF)
  fPsi(i,:) = Lpsi(:,:,i) * psi.' + LdelPsi(:,:,i) * delPsi.';
end

% --- spectrum plot
A= fftshift(fft2([Lpsi(:,:,1),Lpsi(:,:,2)] ));
figure; surf( 1:size(A,2),...
  -floor(size(fPsi,2)/2):(floor(size(fPsi,2)/2)-1) ,...
  20*log10(abs(A)),'EdgeColor','none');
view(135,45)
xlabel('Sphere samples', 'FontSize', 12)
ylabel('DFT spectrum', 'FontSize', 12)
zlabel('Coefficient amplitude [dB]', 'FontSize', 12)
axis tight
% printEPS('','DFTspectrum')

% --- pattern plot
gain = sf_computeGain(fPsi);
gainRef = sf_computeGain(fPsiRef);
sf_plotFFCutPlanes(thetaFF, gain, thetaFF, gainRef, [1 2]);