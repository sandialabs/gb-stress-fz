%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to map a stressed grain boundary into the equivalent
% disorientation, boundary plane and stress fundamental zones.
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
% stress            Stress tensor in Voigt notation
% Optional:  
% 'snap2grid'       Snap points to a boundary plane grid gb_grid, and use 
%                   the symmetries of that point
%
%%% Outputs:
% tests   - Tests structure
%   id              Unique id for each condition [id_dis, id_gb]
%   refgrain        Grain taken as the reference I (1 or 2)
%   id_gbgrid       ID of the grid point closest to each condition (only
%                   for 'snap2grid' option)
%   th_gbgrid       Angle [deg] to the closest grid point (only for
%                   'snap2grid' option)
%   fz              Location in the disorientation FZ
%   fzboundary      Location in the boundary plane FZ
%   fzloadaxis      Location in the load axis FZ
%   symdichromatic  Point group of the boundary plane FZ
%   symgboundary    Point group of the load axis FZ
%   w               Weights (non applicable)
%   euleranglesI    Euler angles of grain I [0,0,0]
%   euleranglesII   Euler angles of grain II (ZXZ)
%   l               Disorientation axis
%   th              Disorientation angle [rad]
%   rv              Rodrigues vector
%   thrangeboundary Angle range of the boundary plane FZ
%   symx,symy,symz  Axes in the symmetry coordinate system
%   phiboundary,thboundary  Spherical coordinates (phi,th) of the boundary planes in the symmetry coordinate system
%   nboundary       Miller indices of the normal to the boundary plane in grain I
%   GBangles        Angles [thB,phiB] of the boundary normal
%   thrangeloadaxis Angle range of the load axis FZ
%   Bx,By,Bz,Bx2    Axes in the load axis coordinate system
%   philoadaxisB,thloadaxisB  Spherical coordinates (phi,th) of the load axis in the symmetry coordinate system
%   nloadaxis       Miller indices of the load axis in grain I coordinate system
%   stressR         Remote stress (Voigt notation) in grain I coordinate system
%

function tests = toFZ_stress(euleranglesI,euleranglesII,GBangles,stress,varargin)
    %%%%% Asserts & initialization   
    assert(size(euleranglesI,1) == size(euleranglesII,1) && ...
           size(euleranglesI ,2) == 3 && ...
           size(euleranglesII,2) == 3,'euleranglesI and euleranglesII must be nx3 matrices');
    assert(size(euleranglesI,1) == size(GBangles,1) && ...
           size(GBangles,2) == 2,'GBangles must be an nx2 matrix');
    assert(size(stress,1) == 1 && ...
           size(stress,2) == 6,'stress must be a 1x6 matrix');

    % Additional options
    snap2grid = false;
    for i = 1:2:length(varargin)                
        propName = varargin{i};
        propValue = varargin{i+1};
        switch lower(propName)
            case 'snap2grid'
                snap2grid = true;
                gb_grid = propValue;
                assert(isstruct(gb_grid),'gb_grid must be a boundary plane structure')
            otherwise
                warning('Unknown property: %s', propName);
        end
    end

    n_in = size(euleranglesI,1);                                    % Number of points
    
    %%%%% Stress symmetries - extracted from the Eigenvalues
    % Determining eigenvalues and eigenvectors
    stressM = trans.V2Nstress(stress);                              % Input stress
    [evectors,evalues] = eig(stressM);                              % Eigenvalues & eigenvectors
    [evalues,isort] = sort(diag(evalues));                          % Sorting
    evectors = evectors(:,isort');
    for j = 1:3                                                     % Normalizing
        evectors(:,j) = evectors(:,j)/norm(evectors(:,j));
    end                                                             % Ensuring right-handedness
    if abs(sum(cross(evectors(:,1),evectors(:,2))-evectors(:,3))) > 10*eps
        evectors(:,1) = -evectors(:,1);
    end
    [evunique,~,ic] = uniquetol(evalues,10*eps);                    % Unique Eigenvalues
    nev = length(evunique);                                         % Number of unique eigenvalues

    % Stress type 
    if nev == 1                                                     % Hydrostatic
        stresstype = 1;        
    elseif nev == 2                                                 % Axial symmetry
        stresstype = 2;
        phimax = pi/2;
        ic = accumarray(ic,1);                                      % Steps to determine the unique eigenvalue
        if ic(1) == 1
            evectors = evectors(:,[2,3,1]);
        else
            evectors = evectors(:,[1,2,3]);
        end
        loadaxis = evectors(:,3)';
    else                                                            % Triaxial
        stresstype = 3;
        evectors = evectors(:,[1,2,3]);
        stress0 = [evalues(1),evalues(2),evalues(3),0,0,0]';
    end      

    %%%%% FZ boundary plane
    [gb,r02gbfz] = toFZ_boundaryplane(euleranglesI,euleranglesII,GBangles);

    % Snapped to disorientation FZ?
    if snap2grid
        [id_gbgrid,th_gbgrid,gb] = snap2grid_boundaryplane(gb,gb_grid);         %%%%% snap2grid_disorientation
    end
    

    %%%%% Finding the FZ boundary plane
    % Preallocating
    stress_FZ = nan(n_in,6);
    if stresstype == 1                                              % Hydrostatic
        fzH = cell(n_in,1);
    elseif stresstype == 2                                          % Axisymmetric
        fzloadaxis = cell(n_in,1);
        [philoadaxis_FZ,thloadaxis_FZ] = deal(nan(n_in,1));
        nloadaxis_FZ = nan(n_in,3);
    else                                                            % Triaxial
        fzstress = cell(n_in,1);
        stressrv_FZ = nan(n_in,3);
    end


    for i = 1:n_in

        %%%%%%%%%%%%%%%%%%%%%%%%%%% HYDROSTATIC %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if stresstype == 1
            fzH{i} = {'H'};

        %%%%%%%%%%%%%%%%%%%%%%%%%% AXISYMMETRIC %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        elseif stresstype == 2
            % Symmetries
            symrot = symmetryrotations(gb.symgboundary{i});             % Proper symmetry rotations
            nsym = length(symrot);                                      % Number of proper symmetry rotations
    
            % Stress - operating on the main eigenvector (loadaxis) and the whole stress 
            loadaxisI = r02gbfz{i} * loadaxis';                         % I - reference grain I coordinate system
            stressI = trans.RVstress(stress',r02gbfz{i});
    
            rfzB2I = gb.Bx(i,:)' * [1,0,0] + ...                        % Rotation matrix (B symmetry -> grain I)
                     gb.By(i,:)' * [0,1,0] + ...
                     gb.Bz(i,:)' * [0,0,1];
            loadaxisfzB = rfzB2I' * loadaxisI;                          % fzBsym coordinate system
            stressfzB = trans.RVstress(stressI,rfzB2I');
    
            % Finding the stress main axis within the FZ
            for j = 1:(2*nsym)                                          % Considering ...
                if j <= nsym                                            % ... proper rotations and ...
                    symrotj = symrot{j};
                else                                                    % ... improper rotations
                    symrotj = -symrot{j-nsym};
                end
                loadaxiss = symrotj*loadaxisfzB;                        % Symmetric loadaxis vectors
        
                % Axis and angle representation
                [~,thloadaxiss,philoadaxiss] = trans.cart2sph(loadaxiss);
                thloadaxiss = mod(thloadaxiss+2*pi,2*pi);
    
                tol = 10*eps; 
                if philoadaxiss >= -tol && philoadaxiss <= pi/2 + tol && ...            % FZ conditions
                   thloadaxiss >= -tol && thloadaxiss <= gb.thrangeloadaxis(i) + tol   % ... considering a tolerance
    
                    % Recovering the rotations angles
                    philoadaxis_FZ(i) = philoadaxiss;
                    thloadaxis_FZ(i) = thloadaxiss;
    
                    nloadaxis_FZ(i,:) = rfzB2I * loadaxiss;              % loadaxis in I coordinate system
                    stress_FZ(i,:) = trans.RVstress(stressfzB,...       % Stress in I coordinate system
                                     rfzB2I*symrotj);
       
                    % Region of the FZ
                    phimax = pi/2;
                    thmax = gb.thrangeloadaxis(i);
                    fzloadaxis{i} = regionS2(philoadaxis_FZ(i),thloadaxis_FZ(i),phimax,thmax,gb.symgboundary{i});
    
                    % Exit loop
                    break
                end
            end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%% TRIAXIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else

            % Symmetries
            Sgb = symmetryrotations(gb.symgboundary{i});                % Grain boundary - proper rotations
            Ss  = symmetryrotations('D2h');                             % Stress - proper rotations
%             Sgb = [Sgb,cellfun(@(x) -x,SgbP,'UniformOutput',false)];   % Proper & improper
%             Ss = [Ss,cellfun(@(x) -x,SsP,'UniformOutput',false)];
            nsymgb = length(Sgb);                                       % Number of symmetry rotations
            nsyms = length(Ss);
            [V,~] = fzSO3(gb.symgboundary{i},'D2h');                    % Rotation fundamental zone (vertices)

            % Remote stress
            stressI = trans.RVstress(stress',r02gbfz{i});               % I - reference grain I coordinate system
    
            rfzB2I = gb.Bx(i,:)' * [1,0,0] + ...                        % Rotation matrix (B symmetry -> grain I)
                     gb.By(i,:)' * [0,1,0] + ...
                     gb.Bz(i,:)' * [0,0,1];
            stressfzB = trans.RVstress(stressI,rfzB2I');                % fzBsym coordinate system

            % Eigenvalues & eigenvectors of the remote stress stressfzB
            [evectorsfzB,evaluesfzB] = eig(trans.V2Nstress(stressfzB)); % Eigenvalues and eigenvectors
            [~,isortfzB] = sort(diag(evaluesfzB));                      % Sorting
            evectorsfzB = evectorsfzB(:,isortfzB');
            for j = 1:3                                                 % Normalizing
                evectorsfzB(:,j) = evectorsfzB(:,j)/norm(evectorsfzB(:,j));
            end    
                                                                        % Ensuring right-handedness
            if abs(sum(cross(evectorsfzB(:,1),evectorsfzB(:,2))-evectorsfzB(:,3))) > 10*eps
                evectorsfzB(:,1) = -evectorsfzB(:,1);
            end

            % Rotation matrix (stressfzB = rstressB*stress0*rstressB^(-1))
            rstressB = evectorsfzB(:,1) * [1,0,0] + ...                % Active rotation matrix
                       evectorsfzB(:,2) * [0,1,0] + ...
                       evectorsfzB(:,3) * [0,0,1];
                    
            % Finding the Rodrigues vector within the FZ
            tol = 10*eps;                                               % Tolerance
            for j = 1:nsymgb                                            % Considering ...
                for k = 1:nsyms
                    rstressBs = Sgb{j}*rstressB*Ss{k};                  % Symmetric rotation
                    
                    if abs(det(rstressBs)-1) <= tol                     % Axis angle representation
                        stressrvs = trans.r2rodrigues(rstressBs);
                    elseif abs(det(rstressBs)+1) <= tol
                        stressrvs = trans.r2rodrigues(-rstressBs);
                    else
                        error('Determinant of rstressB must be 1 or -1')
                    end
                    
                    if inhull(stressrvs,V,[],10*eps)                    % Is the rodrigues vector within FZ?
       
                        % Recovering the FZ rotations and stresses
                        stressrv_FZ(i,:) = stressrvs;                   % Rodrigues vector within the FZ
                        stressFZB = trans.RVstress(stress0,rstressBs);  % Stress in the B coordinate system and within stress FZ
        
                        stress_FZ(i,:) = trans.RVstress(stressFZB,...   % Stress in I coordinate system
                                         rfzB2I);
               
                        % Region of the FZ
                        ge = 0;                                         % Grain exchange symmetry OFF
                        fzstress{i} = regionSO3(gb.symgboundary{i},'D2h',stressrv_FZ(i,:),ge);

                        % Exit loop
                        break
                    end
                end
            end
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

    end
    
    
    % Building ID - [disorientation,gbplane,loadaxis]
    [~,~,id1] = uniquetol([gb.euleranglesI,gb.euleranglesII,gb.rv],100*eps,'ByRows',true);
    id2 = nan(n_in,1);
    for i = 1:max(id1)
        fu = id1 == i;
        [~,~,id2(fu)] = uniquetol(gb.GBangles(fu,:),100*eps,'ByRows',true);
    end
    [~,~,uid12] = uniquetol([id1,id2],100*eps,'ByRows',true);
    id3 = cell2mat(splitapply(@(x){1:numel(x)}, 1:numel(uid12), uid12'))';

    %%%%% Building dis structure
    tests.id = [id1,id2,id3];
    tests.refgrain = gb.refgrain;
    if snap2grid
        tests.id_gbgrid = id_gbgrid;
        tests.th_gbgrid = th_gbgrid;
    end
    tests.fzdichromatic = gb.fzdichromatic;
    tests.fzboundary = gb.fzboundary;
    if stresstype == 1
        tests.fzH = fzH;
    elseif stresstype == 2
        tests.fzloadaxis = fzloadaxis;
    else
        tests.fzstress = fzstress;
    end
    tests.symdichromatic = gb.symdichromatic;
    tests.symgboundary = gb.symgboundary;
    tests.w = ones(n_in,3);
    tests.euleranglesI = gb.euleranglesI;
    tests.euleranglesII = gb.euleranglesII;
    tests.l = gb.l;
    tests.th = gb.th;
    tests.rv = gb.rv;
    tests.thrangeboundary = gb.thrangeboundary;
    tests.symx = gb.symx;
    tests.symy = gb.symy;
    tests.symz = gb.symz;
    tests.symx2 = gb.symx2;
    tests.phiboundary = gb.phiboundary;
    tests.thboundary = gb.thboundary;
    tests.nboundary = gb.nboundary;
    tests.GBangles = gb.GBangles;
    tests.Bx = gb.Bx;
    tests.By = gb.By;
    tests.Bz = gb.Bz;
    tests.Bx2 = gb.Bx2;
    tests.thrangeloadaxis = gb.thrangeloadaxis;
    if stresstype == 2
        tests.philoadaxisB = philoadaxis_FZ;
        tests.thloadaxisB = thloadaxis_FZ;
        tests.nloadaxis = nloadaxis_FZ;
    elseif stresstype == 3
        tests.stressrv = stressrv_FZ;
    end
    tests.stressR = stress_FZ;
end