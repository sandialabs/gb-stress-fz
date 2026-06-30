%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating elastic incompatibility stresses at a stressed grain boundary
% and related damage metrics.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculation of the incompatibility stress is performed in ss_bicrystal.m
% using the model by Richeton et al. (2015).
%
% Coordinate systems (CS):
% _T            Global test
% _B            Boundary (y//GB normal))
% _I            Grain I - crystallographic
% _II           Grain II - crystallographic
%
%%% Inputs
% scomponents   Compliance components [s011,s012,s044]
% fI            Volume fraction of grain I
% stressRT      Remote stress tensor in Voigt notation
% euleranglesI  Euler angles - grain I  [phi1I ,PhiI ,phi2I ]
% euleranglesII Euler angles - grain II [phi1II,PhiII,phi2II]
% GBangles      Grain boundary angles [thetaB,phiB] (polar and azimuthal)
% outputs       Outputs for sI and sII:
%               'mS'  -> mS and stressmag
%               'sp' -> add slip plane (sp.[taup,phi,mN,mS]) description,
%                       only if 'mS' included
%               'mT' -> add table T, only if 'mS' included.
%               'mScosphi' -> add mS*cos(phi) included, where phi is the 
%                             angle between the GB and the slip plane.
%               'strainenergy' -> add strain energy density
%               'all' -> add all of the above
%
%%% Outputs:
% sI, sII       Outputs for grains I and II
%  mS           List of Schmid factors
%  stressmag    Magnitude of the stress max(ps)-min(ps), such that the
%               resolved shear stress trss = mS*stressmag
%  sp           Slip plane description (see Leon-Cazares et al. 2020, Crystals)
%  mT           Table with all 
%  mScosphi     mS*cos(phi); phi is the angle between boundary and slip plane normals
%  strainenergy Strain energy density [J/m^3]
% stressCB      Incompatibility stress at the GB, in B CS
%

function [sI,sII,stressCB] = bicrystal(scomponents,fI,stressRT,euleranglesI,euleranglesII,GBangles,outputs)

    assert(all(size(scomponents)==[1,3]),'scomponents must be a 1x3 vector')
    assert(all(size(fI)==[1,1]),'fI must be a scalar value')
    assert(all(size(stressRT)==[6,1]),'stressRT must be a 6x1 vector')
    assert(all(size(euleranglesI)==[1,3]),'euleranglesI must be a 1x3 vector')
    assert(all(size(euleranglesII)==[1,3]),'euleranglesII must be a 1x3 vector')
    assert(all(size(GBangles)==[1,2]),'GBangles must be a 1x2 vector')
    
    [sI,sII] = deal([]);        % Pre-assigning
    if strcmp(outputs,'all')    % All outputs
        outputs = {'mS','sp','mT','mScosphi','strainenergy'};
    end
    
    %%%%% Inputs
    s011 = scomponents(1);      % Compliance components
    s012 = scomponents(2);
    s044 = scomponents(3);
    
    phi1I  = euleranglesI(1);   % Euler angles - grain I
    PhiI   = euleranglesI(2);
    phi2I  = euleranglesI(3);
    
    phi1II = euleranglesII(1);  % Euler angles - grain II
    PhiII  = euleranglesII(2);
    phi2II = euleranglesII(3);
    
    thetaB = GBangles(1);       % Grain boundary angles (polar and azimuthal)
    phiB = GBangles(2);
    
    %%%%% Defining the CSs
    % Defining the boundary plane in the T CS
    nboundaryT = trans.sph2cart(1,thetaB,phiB);       % Normal to the GB
    if nboundaryT(3) < 1                                % Two vectors perpendicular to nboundary
        nboundaryT_perp1 = [ nboundaryT(2)/sqrt(nboundaryT(1)^2+nboundaryT(2)^2),...
                            -nboundaryT(1)/sqrt(nboundaryT(1)^2+nboundaryT(2)^2),...
                             0]';                        %   _perp1 chosen perpendicular to z
    else
        nboundaryT_perp1 = [1,0,0]';
    end
    nboundaryT_perp2 = cross(nboundaryT_perp1,nboundaryT);
    
    % Rotation matrices: changing vectors from X to Y such that RX2Y * _X = _Y
    rT2I  = trans.rEuler(phi1I ,PhiI ,phi2I );
    rT2II = trans.rEuler(phi1II,PhiII,phi2II);
    rT2B = [1;0;0]*nboundaryT_perp1' + ...   
           [0;1;0]*nboundaryT' + ...
           [0;0;1]*nboundaryT_perp2';
    rB2I  = rT2I *(rT2B)^-1;                              
    rB2II = rT2II*(rT2B)^-1;                              
    
    %%%%% Compliance
    s0 = [s011, s012, s012,    0,    0,    0;...        % Cubic compliance matrix in the crystal CS (I & II), Voigt notation
          s012, s011, s012,    0,    0,    0;...
          s012, s012, s011,    0,    0,    0;...
             0,    0,    0, s044,    0,    0;...
             0,    0,    0,    0, s044,    0;...
             0,    0,    0,    0,    0, s044];
    sIB  = trans.RVcompliance(s0,rB2I ^-1);             % Compliance of each grain in Boundary CS, Voigt notation
    sIIB = trans.RVcompliance(s0,rB2II^-1);
    
    
    %%%%% Incompatibility stress at the GB
    % Tensors in the B CS  
    stressRB  = trans.RVstress(stressRT,rT2B);          % Remote stress in B
    
    % [Richeton 2015, Phil. Mag. 95:1]; stress in each grain at the GB in the B CS
    [stressBI_B,stressBII_B,stressCB] = ss_bicrystal(sIB,sIIB,fI,stressRB,0);
    
    stressBI_I   = trans.RVstress(stressBI_B ,rB2I );   % Total stress in each grain, in grain CS
    stressBII_II = trans.RVstress(stressBII_B,rB2II);
    
    %%%%% Resolved shear stresses
    % [Leon-Cazares 2020, Crystals 10]; Schmid, Escaig and Normal factor, and slip plane descriptors
    if any(strcmp(outputs,'mS'))
        [ sI.mS, sI.mE, sI.stressmag, sI.sp, sI.mT] = ss_slipsystems(stressBI_I  ,12,outputs);
        [sII.mS,sII.mE,sII.stressmag,sII.sp,sII.mT] = ss_slipsystems(stressBII_II,12,outputs);
    end
    
    % Calculating the products mS*cos(phi), where phi is the angle between the
    % boundary plane normal and the slip plane normal.
    if any(strcmp(outputs,'mScosphi'))
        n =  [ 1, 1, 1; 1, 1, 1; 1, 1, 1;                   % Slip plane normals in I & II CS (matching ss_slipsystems nomenclature)
              -1,-1, 1;-1,-1, 1;-1,-1, 1;
               1,-1,-1; 1,-1,-1; 1,-1,-1;
              -1, 1,-1;-1, 1,-1;-1, 1,-1]/sqrt(3);
        nboundaryI  = rB2I *[0;1;0];                        % Boundary plane in I & II CS
        nboundaryII = rB2II*[0;1;0];
        [sI.mScosphi,sII.mScosphi] = deal(nan(12,1));
        for i = 1:12
            sI.mScosphi(i)  =  sI.mS(i)*dot(nboundaryI ,n(i,:));
            sII.mScosphi(i) = sII.mS(i)*dot(nboundaryII,n(i,:));
        end
    end
    
    %%%%% Strain energy density: U = 1/2*sum(si*ei)
    if any(strcmp(outputs,'strainenergy'))
        s = num2cell(stressBI_I);           % Grain I
        [sxx,syy,szz,syz,sxz,sxy] = deal(s{:});
        sI.strainenergy  = 1/2*(s011*(sxx^2+syy^2+szz^2) + ...
                                 2*s012*(sxx*syy+syy*szz+szz*sxx) + ...
                                 s044*(sxy^2+syz^2+sxz^2));
        s = num2cell(stressBII_II);         % Grain II
        [sxx,syy,szz,syz,sxz,sxy] = deal(s{:});
        sII.strainenergy = 1/2*(s011*(sxx^2+syy^2+szz^2) + ...
                                 2*s012*(sxx*syy+syy*szz+szz*sxx) + ...
                                 s044*(sxy^2+syz^2+sxz^2));
    end

end






