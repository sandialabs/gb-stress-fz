%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to map a grain boundary into the equivalent disorientation 
% and boundary plane in the fundamental zones.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%% Inputs
% euleranglesI      List of euler angles (ZXZ) for each grain
% euleranglesII
% GBangles          Angles [th,phi] of the boundary normal
% *Optional inputs:  
% 'snap2grid'       Snap points to a disorientation grid dis_grid, and use 
%                   the symmetries of that point
%
%%% Outputs:
% gb   - Grain boundary structure
%   id              Unique id for each condition [id_dis, id_gb]
%   refgrain        Grain taken as the reference I (1 or 2)
%   id_disgrid      ID of the grid point closest to each condition (only
%                   for 'snap2grid' option.
%   th_disgrid      Angle [deg] to the closest grid point (only for 
%                   'snap2grid' option.
%   fzdichromatic   Location in the disorientation FZ
%   fzboundary      Location in the boundary plane FZ
%   symdichromatic  Point group of the boundary plane FZ
%   symgboundary    Point group of the grain boundary
%   w               Weights (non applicable)
%   euleranglesI    Euler angles of grain I [0,0,0]
%   euleranglesII   Euler angles of grain II (ZXZ)
%   l               Disorientation axis
%   th              Disorientation angle [rad]
%   rv              Rodrigues vector
%   thrangeboundary Angle range of the boundary plane FZ
%   symx,symy,symz,symx2    Axes in the boundary plane symmetry coordinate system
%   phiboundary,thboundary  Spherical coordinates (phi,th) of the boundary planes in the symmetry coordinate system
%   nboundary       Miller indices of the normal to the boundary plane in grain I coordinate system
%   GBangles        Angles [thB,phiB] of the boundary normal
%   Bx,By,Bz,Bx2    Symmetry axes of the grain boundary plane coordinate system
%   thrangeloadaxis Angle range of the load axis FZ
% r02gbfz           Rotation matrix from the original to the boundary plane coordinate system in the FZ
%

function [gb,r02gbfz] = toFZ_boundaryplane(euleranglesI,euleranglesII,GBangles,varargin)
    %%%%% Asserts & initialization   
    assert(size(euleranglesI,1) == size(euleranglesII,1) && ...
           size(euleranglesI ,2) == 3 && ...
           size(euleranglesII,2) == 3,'euleranglesI and euleranglesII must be nx3 matrices');
    assert(size(euleranglesI,1) == size(GBangles,1) && ...
           size(GBangles,2) == 2,'GBangles must be an nx2 matrix');

    % Additional options
    snap2grid = false;
    for i = 1:2:length(varargin)                
        propName = varargin{i};
        propValue = varargin{i+1};
        switch lower(propName)
            case 'snap2grid'
                snap2grid = true;
                dis_grid = propValue;
                assert(isstruct(dis_grid),'dis_grid must be a disorientation structure')
            otherwise
                warning('Unknown property: %s', propName);
        end
    end

    % Preallocating
    n_in = size(euleranglesI,1);                                    % Number of points
    [phiboundary_FZ,thboundary_FZ] = deal(nan(n_in,1));
    GBangles_FZ = nan(n_in,2);
    nboundary_FZ = nan(n_in,3);
    [r02gbfz,fzboundary] = deal(cell(n_in,1));
    
    %%%%% FZ disorientation
    [dis,r02disfz] = toFZ_disorientation(euleranglesI,euleranglesII);
    
    % Snapped to disorientation FZ?
    if snap2grid
        [id_disgrid,th_disgrid,dis] = snap2grid_disorientation(dis,dis_grid);         %%%%% snap2grid_disorientation
    end


    %%%%% Finding the FZ boundary plane
    for i = 1:n_in

        % Symmetries
        symrot = symmetryrotations(dis.symdichromatic{i});          % Proper symmetry rotations
        nsym = length(symrot);                                      % Number of proper symmetry rotations

        % GBnormal
        nboundary0 = trans.sph2cart(1,GBangles(i,1),GBangles(i,2)); % 0 - initial coordinate system
        nboundaryI = r02disfz{i} * nboundary0;                      % I - reference grain I coordinate system

        rfzsym2I = dis.symx(i,:)' * [1,0,0] + ...                   % Rotation matrix (symmetry -> grain I)
                   dis.symy(i,:)' * [0,1,0] + ...
                   dis.symz(i,:)' * [0,0,1];
        nboundaryfzsym = rfzsym2I' * nboundaryI;                    % GB FZsym coordinate system
 
        % Finding the GBnormal within the FZ
        for j = 1:(2*nsym)                                          % Considering ...
            if j <= nsym                                            % ... proper rotations and ...
                symrotj = symrot{j};
            else                                                    % ... improper rotations
                symrotj = -symrot{j-nsym};
            end
            nboundarys = symrotj*nboundaryfzsym;                    % Symmetric nboundary vectors
    
            % Axis and angle representation
            [~,thboundarys,phiboundarys] = trans.cart2sph(nboundarys);
            thboundarys = mod(thboundarys+2*pi,2*pi);

            tol = 10*eps; 
            if phiboundarys >= -tol && phiboundarys <= pi/2 + tol && ...            % FZ conditions
               thboundarys >= -tol && thboundarys <= dis.thrangeboundary(i) + tol   % ... considering a tolerance

                % Recovering the rotations angles
                phiboundary_FZ(i) = phiboundarys;
                thboundary_FZ(i) = thboundarys;
                
                nboundary_FZ(i,:) = rfzsym2I * nboundarys;              % nboundary in I coordinate system
                [~,thI,phiI] = trans.cart2sph(nboundary_FZ(i,:));       % ... and its corresponding GBangles [thB,phiB]
                GBangles_FZ(i,:) = [thI,phiI];

                % Rotation from 0 to gbFZ coordinate system
                r02gbfz{i} = rfzsym2I * symrotj * rfzsym2I' * r02disfz{i};

                % Region of the FZ
                phimax = pi/2;
                thmax = dis.thrangeboundary(i);
                fzboundary{i} = regionS2(phiboundary_FZ(i),thboundary_FZ(i),phimax,thmax,dis.symdichromatic{i});

                % Exit loop
                break
            end
        end
    end

    % Assigning nan to 'Oh' disorientations (no grain boundary)
    iOh = find(strcmp(dis.symdichromatic,'Oh'));
    for i = iOh'
       fzboundary{i} = 'nan'; 
       r02gbfz{i} = r02disfz{i};
    end

    %%% Miller indices of the boundary normal in grain I cs
    symgboundary = cell(size(nboundary_FZ,1),1);            % Symmetry
    thrangeloadaxis = nan(size(nboundary_FZ,1),1);          % th range for the boundary FZ
    [Bx,By,Bz,Bx2] = deal(nan(size(nboundary_FZ,1),3));     % Symmetry axes
    for i = 1:size(nboundary_FZ,1)
        [symgboundary{i},thrangeloadaxis(i),Bx(i,:),By(i,:),Bz(i,:),Bx2(i,:)] = fccFZ.stressaxissymmetry(nboundary_FZ(i,:),...
                                                                    dis.fzdichromatic{i},dis.symdichromatic{i},fzboundary{i},dis.thrangeboundary(i), ...
                                                                    dis.symx(i,:),dis.symy(i,:),dis.symz(i,:),dis.symx2(i,:),true);
    end 
    
    % Building ID - [disorientation,gbplane]
    [~,~,id1] = unique([dis.euleranglesI,dis.euleranglesII,dis.rv],'rows','stable');
    id2 = cell2mat(splitapply(@(x){1:numel(x)}, 1:numel(id1), id1'))';

    %%%%% Building dis structure
    gb.id = [id1,id2];
    gb.refgrain = dis.refgrain;
    if snap2grid
        gb.id_disgrid = id_disgrid;
        gb.th_disgrid = th_disgrid;
    end
    gb.fzdichromatic = dis.fzdichromatic;
    gb.fzboundary = fzboundary;
    gb.symdichromatic = dis.symdichromatic;
    gb.symgboundary = symgboundary;        
    gb.w = ones(n_in,2);
    gb.euleranglesI = dis.euleranglesI;
    gb.euleranglesII = dis.euleranglesII;
    gb.l = dis.l;
    gb.th = dis.th;
    gb.rv = dis.rv;
    gb.symx = dis.symx;
    gb.symy = dis.symy;
    gb.symz = dis.symz;
    gb.symx2 = dis.symx2;
    gb.thrangeboundary = dis.thrangeboundary;
    gb.phiboundary = phiboundary_FZ;
    gb.thboundary = thboundary_FZ;
    gb.nboundary = nboundary_FZ;
    gb.GBangles = GBangles_FZ;
    gb.Bx = Bx;
    gb.By = By;
    gb.Bz = Bz;
    gb.Bx2 = Bx2;
    gb.thrangeloadaxis = thrangeloadaxis;    
end