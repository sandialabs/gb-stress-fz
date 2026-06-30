%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating incompatibility stresses at a stressed grain boundary.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Incompatibility stresses at a bicrystal, following the derivation by: 
% T. Richeton, I. Tiba, S. Berbenni & O. Bouaziz (2015) Analytical 
% expressions of incompatibility stresses at Sigma3⟨111⟩ twin boundaries and
% consequences on single-slip promotion parallel to twin plane, 
% Philosophical Magazine, 95:1, 12-31, DOI: 10.1080/14786435.2014.984787.
%
% Elastic compliances s and plastic strains ep are assumed to be piecewise
% uniform. The B (boundary) coordinate system is used, set such that the
% interface corresponds to the x1x3 plane and the x2 axis extends into
% crystal II (-x2 extends into % crystal I). Properties are invariant along
% the x1x3 plane.
%
%%% Inputs:
% sI,sII        Compliance tensor in the Boundary coordinate system
% fI            Volume fraction of grain I
% stressR       Remote stress (column vector in Voigt notation)
% strainP       Plastic strain (either 0 or column vectors in Voigt notation {strainPI,strainPII})
%
%%% Outputs:
% stressI,stressII  Stress tensors in each grain
% stressIncS        Incompatibility stress tensor
%

function [stressI,stressII,stressIncS] = ss_bicrystal(sI,sII,fI,stressR,strainP)

% Elastic and plastic parameters
fII = 1 - fI;                                                               % Volume fraction of grain II
                                                                            % Plastic strain allocation
if strainP == 0                                                             %   Zero plastic strain
    [epI,epII] = deal([0,0,0,0,0,0]');                                                      
else                                                                        %   Plastic strain in each grain
    epI = strainP{1};
    epII = strainP{2};
end

% Calculation of the incompatibility stress tensor
sw = fII*sI + fI*sII;                                                       % Mean elastic compliance tensor
d = sw(1,1)*sw(3,5)^2 + sw(3,3)*sw(1,5)^2 + sw(5,5)*sw(1,3)^2 ...
    - sw(1,1)*sw(3,3)*sw(5,5) - 2*sw(1,3)*sw(1,5)*sw(3,5);
G = zeros(6,6);                                                             % Compliance-related tensor
G(1,1) = (sw(3,3)*sw(5,5) - sw(3,5)^2)/d;
G(3,3) = (sw(1,1)*sw(5,5) - sw(1,5)^2)/d;
G(5,5) = (sw(1,1)*sw(3,3) - sw(1,3)^2)/d;
[G(1,3),G(3,1)] = deal((sw(1,5)*sw(3,5) - sw(1,3)*sw(5,5))/d);
[G(1,5),G(5,1)] = deal((sw(1,3)*sw(3,5) - sw(1,5)*sw(3,3))/d);
[G(3,5),G(5,3)] = deal((sw(1,3)*sw(1,5) - sw(3,5)*sw(1,1))/d);

sJ = sII - sI;                                                              % Compliance jump at the interface
epJ = epII - epI;                                                           % Plastic strain jump at the interface
eSJ = sJ * stressR + epJ;                                                   % Total strain jump at the interface

stressIncS = G * eSJ;                                                       % Incompatibility stress tensor

stressI  = stressR - fII*stressIncS;                                        % Stresses in each grain
stressII = stressR + fI *stressIncS;


end
