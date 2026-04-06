%%% low-rank approximation of the near field to far field operator for a
%%% bounding sphere with compact polynomial functions
clear; clc; close all;
addpath('../../matlabLib');

nbrElems_x = 9; % nbr of point sources for a linear array
Nbr = 11:4:31; % nbr of coefficients to retain in the DFT-truncation
Shape = 1; % 1=1stOrder, 2=2ndOrder, 3=3rdOrder, 4=PWS
if Shape == 1
  Nbr = 11:3:201;
elseif Shape == 2
  Nbr = 11:2:201;
elseif Shape == 3
  Nbr = 10:3:201;
elseif Shape == 4
  Nbr = 11:3:201;
end

arrayPos = buildArray(1, nbrElems_x, .5, 1, .5);
[xMin, xMax, yMin, yMax, zMin, zMax, xPts, yPts, zPts] = ...
  getBoxDim(1, arrayPos, .5, .1, .5 );
[boxPos, boxN, dS, mSize] = ...
  buildBox([1 1 1 1 1 1], xMin, xMax, yMin, yMax, zMin, zMax,...
  xPts, yPts, zPts, 1, 0, 0);
[Rmag, NdotRV] = getBoxVectors(arrayPos, boxPos, boxN);
excitPhasor = sf_excitations(1, arrayPos, 0, 0);
[psi, delPsi] = sf_nfSolver(1, excitPhasor, Rmag, NdotRV);
psi = psi.';
delPsi = delPsi.';

L2Error = zeros(1,length(Nbr));
L2ErrorRef = zeros(1,length(Nbr));
Vects = zeros(1,length(Nbr));
for alpha = 1:length(Nbr)
    
  NbrOfVectors = Nbr(length(Nbr)+1-alpha);
  VectorSelection = (1:NbrOfVectors-1);
  Inc = floor(1100/NbrOfVectors);
  Indices = [1,floor(VectorSelection*Inc)];
  TotalNbrOfVectors = Indices(NbrOfVectors);
  TotNbrVects(alpha) = TotalNbrOfVectors;
  avgVects = mean(TotNbrVects);
  fprintf('--> nbr vectors = %2.4g, avg nbr vectors = %2.4g\n', ...
    TotalNbrOfVectors, avgVects);

  thetaFF = linspace(-pi/2,pi/2,TotalNbrOfVectors);
  phiFF = 0;

  [A, B] = sf_nf2ffOperator(1, thetaFF, phiFF, ...
    boxPos, boxN, dS);
  fPsiRef = sf_directffSolver(1, thetaFF, phiFF, ...
    excitPhasor, arrayPos);
  Size = size(A); 
  TotalNbrOfVectors = length(thetaFF);
  
  % ---- build Lq
  BasisA = A(Indices(:),:);
  BasisB = B(Indices(:),:);
  
  %%---- ShapeFunctions -> building Beta Matrix
  if(Shape == 1) %%% 1st Order
    Beta = zeros(length(thetaFF),NbrOfVectors);
    Range = Indices(1):Indices(2);
    Beta(Range,1) = (Range(length(Range)) - (Range))/...
      (Range(length(Range))- Range(1));
    for i=2:NbrOfVectors-1
      Range = Indices(i-1):Indices(i);
      Beta(Range,i) = ( Range - Range(1))/ ...
        (Range(length(Range)) - Range(1));
      Range = Indices(i):Indices(i+1);
      Beta(Range,i) = (Range(length(Range)) - (Range))./...
        (Range(length(Range)) - Range(1));
    end
    Range = Indices(NbrOfVectors-1):Indices(NbrOfVectors);
    Beta(Range, NbrOfVectors) = ( Range - Range(1))/...
      (Range(length(Range)) - Range(1));
  elseif(Shape == 2) %%% 2nd Order
    Beta = zeros(length(thetaFF),NbrOfVectors);
    for i=1:2:NbrOfVectors-2
      Range = Indices(i):Indices(i+2);
      Xrev =(Range(length(Range)) - (Range))/ ...
        (Range(length(Range))- Range(1));    
      Xfor = ( Range - Range(1))/ (Range(length(Range)) - Range(1));
      Beta(Range,i) = Xrev .* (2*Xrev - 1);
      Beta(Range,i+1) = 4*Xrev.*Xfor;
      Beta(Range,i+2) =  Xfor .* (2*Xfor - 1);
    end
  elseif(Shape == 3) %%% 3rd Order
    Beta = zeros(length(thetaFF),NbrOfVectors);
    for i=1:3:NbrOfVectors-3
      Range = Indices(i):Indices(i+3);
      Xfor = ( Range - Range(1))/ (Range(length(Range)) - Range(1));
      Beta(Range,i) = (1-Xfor).* (2-3*Xfor) .* (1-3*Xfor) ./2;
      Beta(Range,i+1) = 9* Xfor .* (1-Xfor) .* (2-3*Xfor) ./2;
      Beta(Range,i+2) = 9* Xfor .* (1-Xfor) .* (3*Xfor-1) ./2;
      Beta(Range,i+3) =  Xfor .* (2-3*Xfor) .* (1-3*Xfor)/2;
    end
  elseif(Shape == 4)
    Beta = zeros(length(thetaFF),NbrOfVectors);
    k=.1*pi;
    for i=1:NbrOfVectors-1
      Range = Indices(i):Indices(i+1);
      Xrev = -(Range(length(Range)) - (Range))/ (Range(length(Range))- Range(1));    
      Xfor = -( Range - Range(1))/ (Range(length(Range)) - Range(1));
      Beta(Range,i+1) = sin(k*Xfor)./sin(k*(Xrev(1)));
      Beta(Range,i) = sin(k*(Xrev))./sin(k*(Xfor(length(Xrev))));
    end
  end
  
  %%---- Compute patterns
  if ~any(Beta(:,size(Beta,2)).')
    disp('Error in the ')
    return;
  end
  approxPattern = Beta*(BasisA* psi) +  Beta*(BasisB * delPsi);
  pattern =  A * psi +  B * delPsi;

  fprintf('(%d)\n', NbrOfVectors);
  Error = getL2error(approxPattern, pattern);
  RefError = getL2error(approxPattern, fPsiRef.');

  L2Error(alpha) = Error;
  L2ErrorRef(alpha) = RefError;
  Vects(alpha) = NbrOfVectors;
end
%% Plot L2 Error on Pattern
figure;
loglog(Vects,L2Error, '-*b', Vects,L2ErrorRef, '-xk', 'LineWidth', 1.5,...
  'MarkerSize', 7);
xlabel('No of samples Q', 'FontSize', 12)
ylabel('Relative Error', 'FontSize', 12)
legend('N2F', 'Direct');
axis tight;
%printEPS('',['errorShape',num2str(Shape)]);
