% Plots the infinitesimal currents that compose an array
%
% vf_plotArrayGeom(arrayPos, lambda, J, M)
%
% IN: lambda = wavelength
%     arrayPos = array elements locations
%     J = electric currents vectors in Cartesian components
%     M = magnetic currents vectors in Cartesian components
%
% Laurent Ntibarikure
function vf_plotArrayGeom(lambda, arrayPos, J, M)

Jv = real(J./mean(sqrt(sum(J.^2,1))))*lambda/10;
Mv = real(M./mean(sqrt(sum(M.^2,1))))*lambda/10;

figure(gcf);
hidden off
hidden on
hold on;
xlabel('x');
ylabel('y');
zlabel('z');
axis tight;
axis equal;
view(45,45);
for i=1:size(arrayPos,2)
  quiver3(arrayPos(1,i), arrayPos(2,i), arrayPos(3,i), Jv(1,i), Jv(2,i), Jv(3,i), ...
    1, 'b', 'LineWidth',2);
  hold on;
  quiver3(arrayPos(1,i), arrayPos(2,i), arrayPos(3,i), Mv(1,i), Mv(2,i), Mv(3,i), ...
    1, 'r', 'LineWidth',2);
end