%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to calculate the Schmid stresses on the a/2<110>{111} slip systems
% of an fcc crystal.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Adapted from 'Stress orientation maps', 
% https://github.com/ferleoncazares/Stress-orientation-maps-fcc
%
% Slip planes analysed in the next order: (111), (-1-11), (1-1-1), (-11-1)
% Order of partials to ensure ISF formation and use of the same line vector
%
% 3 coordinate systems considered (CSsample -R1-> CScrystal -R2-> CSss): 
%   - CSsample: CS so that a uniaxial force is applied along [1,0,0]
%   - CScrystal: Miller indices
%   - CSss: CS aligned for each slip system so that
%           b = [1,0,0], cross(b,n) = [0,1,0], n = [0,0,1]
%           (b = slip direction, n = slip plane normal)
%
%%% Inputs:
% Scrystal      Stress tensor in Voigt notation
% m             12 or 24 slip systems
% outputs       'sp' -> add slip plane (sp.[taup,phi,mN,mS]) description
%               'mT' -> add table T
%
%%% Outputs:
% mS            Schmid factors for the a/2<110> slip systems
% mE            Escaig factors for the a/2<110> slip systems%
% sp            Slip plane description
% - taup        Shear stresses on each (111) plane
% - phi         Angle [-30deg,30deg] of the shear stress
%               phi = -30 -> pointing towards a leading partial
%               phi = 0   -> pointing towards a perfect dislocation
%               phi = 30  -> pointing towards a trailing partial
%  - mN         Normal factor for each (111) plane
%  - mS         Highest Schmid factor for each (111) plane
% mT            Table with the results
%  - id         Slip system identifier
%  - n          Slip normal
%  - b          Slip direction (perfect a/2<110>{111} dislocations)
%  - mS         Schmid factor
%  - mE         Escaig factor
%  - mN         Normal factor
%  - bp1        Slip direction of leading partials (a/6<112>{111} dislocations)
%  - bp2        Slip direction of trailing partials (a/6<112>{111} dislocations)
%  - mSp1       Schmid factor of bp1 leading partial
%  - mSp2       Schmid factor of bp2 trailing partial
%

function [mS,mE,Smag,sp,mT] = ss_slipsystems(Scrystal,m,outputs)

assert(((m == 12) || (m == 24)),'m must be 12 or 24')
assert(all(size(Scrystal) == [6,1]),'Scrystal must be a 6-by-1 vector')
sp = [];
mT = [];

Scrystal = trans.V2Nstress(Scrystal);   % Voigt to matrix notations

%%%%%%%%%%%%% Calculation
%%% fcc
% Slip system ID
id = (1:m)';

% Glide Planes
n =  [ 1, 1, 1;     % d
       1, 1, 1;
       1, 1, 1;
      -1,-1, 1;     % c
      -1,-1, 1;
      -1,-1, 1;
       1,-1,-1;     % b
       1,-1,-1;
       1,-1,-1;
      -1, 1,-1;     % a
      -1, 1,-1;
      -1, 1,-1];
% Glide Directions
b =  [-1, 0, 1;     % CB
       1,-1, 0;     % BA
       0, 1,-1;     % AC
       1, 0, 1;     % DA
      -1, 1, 0;     % AB
       0,-1,-1;     % BD        
       1, 1, 0;     % DC
       0,-1, 1;     % CA
      -1, 0,-1;     % AD      
       0, 1, 1;     % DB
       1, 0,-1;     % BC
      -1,-1, 0];    % CD

% Slip direction of the partials: 1 = leading, 2 = trailing
% This analysis considers the partials 1=right and 2=left, 
% with mS>0 inducing a force in the dislocation towards the right.
bp1i= [-2, 1, 1;    % dB
        1,-2, 1;    % dA
        1, 1,-2;    % dC
        2,-1, 1;    % gA
       -1, 2, 1;    % gB
       -1,-1,-2;    % gD        
        1, 2,-1;    % bC
        1,-1, 2;    % bA
       -2,-1,-1;    % bD         
       -1, 1, 2;    % aB
        2, 1,-1;    % aC
       -1,-2,-1];   % aD
bp2i= [-1,-1, 2;    % Cd
        2,-1,-1;    % Bd
       -1, 2,-1;    % Ad
        1, 1, 2;    % Dg
       -2, 1,-1;    % Ag
        1,-2,-1;    % Bg         
        2, 1, 1;    % Db
       -1,-2, 1;    % Cb
       -1, 1,-2;    % Ab
        1, 2, 1;    % Da
        1,-1,-2;    % Ba
       -2,-1, 1];   % Ca   

if m == 24
    n = [n;-n];         % By artifitially inverting nsp rather than nsd -> (mS13 = -mS1 & mE13 = mE1) 
    b = [b;b];
    bp1 = [bp1i;bp2i];
    bp2 = [bp2i;bp1i];
else
    bp1 = bp1i;
    bp2 = bp2i;
end

% Determining Scrystal
E = sort(eig(Scrystal));
Smag = (max(E)-min(E));
Scrystaln = Scrystal/Smag;

% Schmid and Escaig factors
mS = zeros(m,1);
mE = zeros(m,1);
mN = zeros(m,1);
mSp1 = zeros(m,1);
mSp2 = zeros(m,1);
for i = 1:m  
    % Rotation matrix: CScrystal -R2-> CSss
    w = cross(n(i,:),b(i,:));
    R2 = [1;0;0]*b(i,:)/sqrt(2) + [0;1;0]*w/sqrt(6) + [0;0;1]*n(i,:)/sqrt(3);
    Sss = R2*Scrystaln*R2';                                             % Stress state in CSss
    mS(i) = Sss(1,3);                                                   % Schmid factor
    mE(i) = Sss(2,3);                                                   % Escaig factor
    if any(strcmp(outputs,'sp')) || any(strcmp(outputs,'mST'))
        if mS(i) < 0                                                        % Reassigning bp1 and bp2 vectors in case mS<0 -> dislocations i and i+12 will now have the same mSp1
            [bp1(i,:),bp2(i,:)] = deal(-bp2(i,:),-bp1(i,:));
        end
        mSp1(i) = dot(Scrystaln*n(i,:)'/sqrt(3),bp1(i,:)'/sqrt(6));         % Schmid factor for each partial dislocation: mSp1 + mSp2 = mS
        mSp2(i) = dot(Scrystaln*n(i,:)'/sqrt(3),bp2(i,:)'/sqrt(6));
        mN(i) = Sss(3,3);                                                   % Stress normal to the slip plane
    end
end
% mS(abs(mS)<1e-15) = 0;          % To ignore small numerical deviations from 0
% mE(abs(mE)<1e-15) = 0;
% mS2(i) = (sum(F'.*n(i,:),2)/sqrt(3) .* sum(F'.*b(i,:),2)/sqrt(2))./sum(F'.^2);    % Gives the same result, directly from the calculations of the individual cosines
% mS3(i) = dot(Scrystal*n(i,:)'/sqrt(3),b(i,:)'/sqrt(2));                           % ",                              from the defintion of the Schmid factor m=Sn.b 


% Slip plane [taup,phi,mN,mS] description
if any(strcmp(outputs,'sp'))           
    [sp.taup,sp.phi,sp.mN] = deal(zeros(4,1));
    np = reshape(1:12,3,4)';        % Indicates which slip systems belong to each slip plane
    for i=1:4
        sp.taup(i) = sqrt(mS(np(i,1)).^2+mE(np(i,1)).^2);
        [~,I] = max(abs(mS(np(i,:))));
        sp.phi(i) = atan2d(mE(np(i,I)),abs(mS(np(i,I))));
        sp.mN(i) = mN(np(i,1));
    end  
    sp.phi(abs(sp.phi)>30.0000001) = NaN;
    sp.mS = sp.taup .* cosd(sp.phi);
end

% Summary table
if any(strcmp(outputs,'mT'))           
    mT = table(id,n,b,mS,mE,mN,bp1,bp2,mSp1,mSp2);              % Table with all the results 
end

end