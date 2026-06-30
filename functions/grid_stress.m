%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to create bicrystal-stress grids of Oh grain boundaries.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Inputs:
% gb                Grain boundary structure:
%   id              Unique id for each condition [id_dis, id_gb]
%   fzdichromatic   Location in the disorientation FZ
%   fzboundary      Location in the boundary plane FZ
%   symdichromatic  Point group of the bicrystal
%   symgboundary    Point group of the grain boundary
%   w               Weights (times each point repeats over the whole orientation space)
%   euleranglesI    Euler angles of grain I [0,0,0]
%   euleranglesII   Euler angles of grain II (ZXZ)
%   l               Disorientation axis
%   th              Disorientation angle [rad]
%   rv              Rodrigues vector
%   symx,symy,symz,symx2    Axes in the symmetry coordinate system
%   thrangeboundary         Angle range of the boundary plane FZ
%   phiboundary,thboundary  Spherical coordinates (phi,th) of the boundary planes in the symmetry coordinate system
%   nboundary       Miller indices of the normal to the boundary plane in grain I
%   GBangles        Angles [thB,phiB] of the boundary normal
%   Bx,By,Bz,Bx2    Symmetry axes of the grain boundary plane coordinate system
%   thrangeloadaxis Angle range of the load axis FZ
% type              Stress mode: 'hydrostatic', 'uniaxial', 'biaxial', 'triaxial', 'shear' or 'tensor'
% stress            Stress values:
%                       type        stress
%                   'hydrostatic'   [sh]
%                   'uniaxial'      [s3]
%                   'biaxial'       [s1,s3] (s2 = s1)
%                   'triaxial'      [s1,s2,s3]
%                   'shear'         [s13]
%                   'tensor'        [s11,s22,s33,s23,s13,s12]
% steps             Number of steps for each degree of freedom
%
% *Optional inputs:
% 'symmetry'        Consider crystal symmetry to create the grid?
%                   false -> no; 2*pi range
%                   true (default) -> yes; only fundamental zone
%                   2 -> yes; only fundamental zone but do not delete
%                        edge points that appear multiple times
% 'txt'             false (default) or true. Write arrays to file to avoid 
%                   preallocation, slower but less RAM.
%
% Outputs:
% tests - Set of tests (grain boundary + load)
%   id              [Disorientation id, grain boundary id, load id]
%   []              All inputs from gb
%   Bx,By,Bz,Bx2    Symmetry axes of the grain boundary plane coordinate system
%   stressR         Remote stress (Voigt notation) in grain I coordinate system
%   Depending on the stresstype:
%   - Hydrostatic stress:
%     fzH           Marker to indicate that this is a hydrostatic stress
%   - Axisymmetric stress:
%     fzloadaxis    Location in the axisymmetric stress FZ
%     philoadaxisB, thloadaxisB   Spherical coordinates (phi,th) of the load axis in the B coordinate system
%     nloadaxis      Load direction in grain I coordinate system
%   - Triaxial stress:
%     fzstress      Location in the triaxial stress FZ
%     stressrv      Stress Rodrigues vector in the B coordinate system
%

function tests = grid_stress(varargin)

    [gb,type,stress,steps] = varargin{1:4};
    % Defaults
    sym = 1;
    txt = false;
    % Additional options
    if length(varargin) > 4
        for i = 5:2:length(varargin)      
            optionName = varargin{i};
            optionValue = varargin{i+1};         

            switch lower(optionName)    % Symmetry
                case 'symmetry'
                    if optionValue == false
                        sym = 0;
                    elseif optionValue == true
                        sym = 1;
                    elseif strcmp(optionValue,'includeedges')
                        sym = 2;
                    else
                        error("symmetry must be false, true (default) or 'includeedges'")  % sym assert
                    end

                case 'txt'
                    if optionValue == false
                        txt = false;
                    elseif optionValue == true
                        txt = true;
                    else
                        error("txt must be false (default) or true")  % txt assert
                    end

                otherwise
                    error('processOptions:UnrecognizedOption', ...
                          'Unrecognized option: ''%s''', optionName);
            end
        end
    end

    % Asserts
    assert(strcmp(type,'hydrostatic') || strcmp(type,'uniaxial') || strcmp(type,'biaxial') || ...
           strcmp(type,'triaxial') || strcmp(type,'shear') || strcmp(type,'tensor'), ...
           "type must be 'hydrostatic', 'uniaxial', 'biaxial', 'triaxial, 'shear' or 'tensor'")
    assert(isnumeric(stress),'Stress must be numeric')
    if strcmp(type,'hydrostatic')
        assert(length(stress)==1,"if type='hydrostatic', then stress must be a scalar value")
    elseif strcmp(type,'uniaxial')
        assert(length(stress)==1,"if type='uniaxial', then stress must be a scalar value")
    elseif strcmp(type,'biaxial')
        assert(all(size(stress)==[2,1]),"if type='biaxial', then stress must be a 2x1 vector")
    elseif strcmp(type,'triaxial')
        assert(all(size(stress)==[3,1]),"if type='triaxial', then stress must be a 3x1 vector")
    elseif strcmp(type,'shear')
        assert(length(stress)==1,"if type='shear', then stress must be a scalar value")
    elseif strcmp(type,'tensor')
        assert(all(size(stress)==[6,1]),"if type='tensor', then stress must be a 6x1 vector")
    end
    assert(isstruct(gb) && isfield(gb,'l') && isfield(gb,'fzdichromatic') && isfield(gb,'GBangles'),'gb must be a structure with specific fields (see grid_boundaryplane.m)')
    assert((length(steps)==1) && all(floor(steps)==steps),'steps must be a single integer')
    
    %%% Number of grain boundaries and grid steps
    gbsize = size(gb.l,1);

    %%% Stress tensor symmetries
    % Stress type
    if strcmp(type,'hydrostatic')
        stressL = [1,1,1,0,0,0]'*stress;
        if steps > 0
            steps = deal(0);
            warning('hydrostatic stress -> steps set to 0')
        end
    elseif strcmp(type,'uniaxial')
        stressL = [0,0,1,0,0,0]'*stress;                    % s33
    elseif strcmp(type,'biaxial')
        stressL = [stress(2),stress(2),stress(1),0,0,0]';   % [s11,s11,s33]
    elseif strcmp(type,'triaxial')
        stressL = [stress(1),stress(2),stress(3),0,0,0]';   % s11, s22, s33
    elseif strcmp(type,'shear')
        stressL = [0,0,0,0,1,0]'*stress;                    % s13
    elseif strcmp(type,'tensor')
        stressL = stress;
    end

    % Stress symmetries - extracted from the Eigenvalues
    stressM = trans.V2Nstress(stressL);         % Stress matrix
    [evectors,evalues] = eig(stressM);          % Eigenvalues and eigenvectors
    [evalues,isort] = sort(diag(evalues));      % Sorting
    evectors = evectors(:,isort');
    for i = 1:3                                 % Normalizing
        evectors(:,i) = evectors(:,i)/norm(evectors(:,i));
    end                                         % Ensuring right-handedness
    if abs(sum(cross(evectors(:,1),evectors(:,2))-evectors(:,3))) > 10*eps
        evectors(:,1) = -evectors(:,1);
    end
    [evunique,~,ic] = uniquetol(evalues,4*eps); % Unique Eigenvalues

    % Rearranging stress
    if length(evunique) == 1                    % Hydrostatic
        stresstype = 1;                       
    elseif length(evunique) == 2                % Axial symmetry
        ic = accumarray(ic,1);
        if ic(1) == 1
            evectors = evectors(:,[2,3,1]);
        else
            evectors = evectors(:,[1,2,3]);
        end
        phimax = pi/2;                          % Range to rotate angle phi (thrangeloadaxis previously defined)
        stresstype = 2;
    else                                        % Triaxial
        evectors = evectors(:,[1,2,3]);
        stresstype = 3;
    end

    % Rotated stress tensor - main symmetry axis parallel to z
    rlstress = evectors(:,1) * [1,0,0] + ...   % Passive rotation matrix
               evectors(:,2) * [0,1,0] + ...
               evectors(:,3) * [0,0,1];
    
    stress0 = trans.RVstress(stressL,rlstress');

    % Load axis symmetry range & symmetry considered
    if sym >= 1                                             % thrangeboundary 
        thrangeloadaxis = gb.thrangeloadaxis;               % Using grain boundary symmetry
        symgboundary = gb.symgboundary;                     % Grain boundary symmetry
    else                                                    
        thrangeloadaxis = zeros(size(gb.thrangeloadaxis)) + 2*pi;   % Whole sircle
        symgboundary = repmat({'Ci'},size(gb.symgboundary,1),1);    % No symmetry considered
    end

    %%% Load grid
    Nloads = nan(gbsize,1);
    sizemax = gbsize*4*max(steps,1)^3;
    if ~txt                     % Normal arrays
        if stresstype == 1              
            stressR = nan(sizemax,6);
            fzH = cell(sizemax,1);
            w = nan(sizemax,1);
        elseif stresstype == 2
            [w,philoadaxisB,thloadaxisB] = deal(nan(sizemax,1));
            fzloadaxis = cell(sizemax,1);
            nloadaxis = nan(sizemax,3);
            stressR = nan(sizemax,6);
        else
            w = nan(sizemax,1);
            fzstress = cell(sizemax,1);
            stressrv = nan(sizemax,3);
            stressR = nan(sizemax,6);
        end
    else                        % Tall arrays - no preallocation
        fprec = 15;                 % Precision (decimal points) for floats
        fw = ['%.',num2str(fprec),'f'];
        f.id = fopen('temp/id.txt','w+');
        f.w = fopen('temp/w.txt','w+');
        if stresstype == 1
            f.stressR = fopen('temp/stressR.txt','w+');
            f.fzH = fopen('temp/fzH.txt','w+');
        elseif stresstype == 2
            f.philoadaxisB = fopen('temp/philoadaxisB.txt','w+');
            f.thloadaxisB = fopen('temp/thloadaxisB.txt','w+');
            f.fzloadaxis = fopen('temp/fzloadaxis.txt','w+');
            f.nloadaxis = fopen('temp/nloadaxis.txt','w+');
            f.stressR = fopen('temp/stressR.txt','w+');
        else
            f.fzstress = fopen('temp/fzstress.txt','w+');
            f.stressrv = fopen('temp/stressrv.txt','w+');
            f.stressR = fopen('temp/stressR.txt','w+');
        end
    end
    ctot = 0;   % Counter     

    % Loop for each GB
    for i = 1:gbsize        

        %%%%%%%%%%%%%%%%%%%%%%%%%%% HYDROSTATIC %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if stresstype == 1
            stressRi = stress0';
            fzHi = {'H'};
            wi = 1;

        %%%%%%%%%%%%%%%%%%%%%%%%%% AXISYMMETRIC %%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif stresstype == 2
            %%% Grid
            thmax = thrangeloadaxis(i);                 % Range of the load axis in S^2
            [phii,thi] = gridS2(phimax,thmax,steps,[]); % Grid over S^2
             
            if all(thi==0), phii = 0; thi = 0; end
            if isnan(thi), thi = 0; end
            
            %%% Nature of the point in the boundary load axis fundamental zone
            fzloadaxisi = cell(length(phii),1);
            wi = nan(length(phii),1);
            for j = 1:length(phii)
                [fzloadaxisi{j},wi(j)] = regionS2(phii(j),thi(j),phimax,thmax,gb.symgboundary{i});
            end
    
            %%% Load axis and stress - B (boundary) coordinate system
            nloadaxisBi = [sin(phii).*cos(thi),sin(phii).*sin(thi),cos(phii)];  % In B coordinate system
                    
            nloadaxisi = nan(length(phii),3);
            stressRi = nan(length(phii),6); 
    
            rB2I = gb.Bx(i,:)' * [1,0,0] + ...                          % Rotation matrix (B -> grain I)
                   gb.By(i,:)' * [0,1,0] + ...
                   gb.Bz(i,:)' * [0,0,1];
    
            for j = 1:length(phii)   
                % Geometry in B coordinate system
                if ~all(nloadaxisBi(j,:) == [0,0,1])
                    load_axis = trans.normalize(cross([0,0,1],nloadaxisBi(j,:)));   % Stress rotation axis 
                    load_angle = acos(dot([0,0,1],nloadaxisBi(j,:)));               % Stress rotation angle
                else
                    load_axis = [0,0,0];
                    load_angle = 0;
                end
                rstressB = trans.raxisangle(load_axis,-load_angle);     % Rotation matrix
                stressRB = trans.RVstress(stress0,rstressB);            % Stress tensor
    
    
                % Geometry in I coordinate system
                nloadaxisi(j,:) = (rB2I * nloadaxisBi(j,:)')';          % Load direction
                stressRi(j,:) = trans.RVstress(stressRB,rB2I);          % Stress tensor
            end
    
            % Deleting repeated points from FZs without a symmetry plane 
            % perpendicular to the main symmetry axis (z) have 
            if sym == 1
                del = zeros(size(fzloadaxisi));
                % C2h symmetry, i.e. viewed along the C2 rotation axis
                % points along theta=180 deg (OB & B)
                % -> equivalent point along theta=0 deg (OA & A)
                del = del | strcmp(gb.symgboundary{i},'C2h') & ...
                            (strcmp(fzloadaxisi,'OB') | strcmp(fzloadaxisi,'B'));
        
                % C2h* symmetry, i.e. C2h viewed along the mirror plane
                % points B and AB for theta>90 deg
                % -> equivalent point along theta=0 deg (OA & A)
                del = del | strcmp(gb.symgboundary{i},'C2h*') & ...
                            ((strcmp(fzloadaxisi,'B')) | ...
                             (strcmp(fzloadaxisi,'AB') & thi>pi/2+8*eps));
        
                % D3d symmetry (OE points), points B and AB for theta>30 deg
                % -> equivalent point at AB for theta < 30 deg
                del = del | strcmp(gb.symgboundary{i},'D3d') & ...
                            ((strcmp(fzloadaxisi,'B')) | ...
                             (strcmp(fzloadaxisi,'AB') & thi>pi/6+8*eps));
                
                % Deleting the points
                phii(del,:) = [];
                thi(del,:) = [];
                fzloadaxisi(del,:) = [];
                wi(del,:) = [];
                nloadaxisi(del,:) = [];
                stressRi(del,:) = [];
            end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%% TRIAXIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            %%% Grid
            ge = 0;                                         % Grain exchange symmetry OFF
            reg = 1;                                        % Output region ON
            psgb = symgboundary{i};                         % Point group of the grain boundary    

            [stressrvi,fzstressi,wi] = gridSO3(psgb,'D2h',steps,ge,reg);   % Grid over SO(3)
                    
            %%% Stress - B (boundary) coordinate system
            stressRi = nan(size(stressrvi,1),6); 
    
            rB2I = gb.Bx(i,:)' * [1,0,0] + ...                          % Rotation matrix (B -> grain I)
                   gb.By(i,:)' * [0,1,0] + ...
                   gb.Bz(i,:)' * [0,0,1];
    
            for j = 1:size(stressrvi,1)   
                % Geometry in B coordinate system
                rstressB = trans.rodrigues2r(stressrvi(j,:));           % Rotation matrix
                stressRB = trans.RVstress(stress0,rstressB);            % Stress tensor
    
                % Geometry in I coordinate system
                stressRi(j,:) = trans.RVstress(stressRB,rB2I);          % Stress tensor
            end
    
            % Deleting repeated points from FZs without a symmetry plane 
            % perpendicular to the main symmetry axis (z) have 
            if sym == 1
               % TODO
            else    % includeedges -> the weights are no longer accurate
                wi = nan(size(wi));
            end               
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Conglomerating and saving
        c = size(stressRi,1);                      % Counter
        Nloads(i) = c;                              % Number of stresses per grain boundary
        if ~txt                                     % Saving values - normal arrays
            if stresstype == 1
                stressR(ctot+1:ctot+c,:) = stressRi;
                fzH(ctot+1:ctot+c) = fzHi;
                w(ctot+1:ctot+c) = wi;
            elseif stresstype == 2
                philoadaxisB(ctot+1:ctot+c) = phii;
                thloadaxisB(ctot+1:ctot+c) = thi;
                fzloadaxis(ctot+1:ctot+c) = fzloadaxisi;
                w(ctot+1:ctot+c) = wi;
                nloadaxis(ctot+1:ctot+c,:) = nloadaxisi;
                stressR(ctot+1:ctot+c,:) = stressRi; 
            else
                fzstress(ctot+1:ctot+c) = fzstressi;
                w(ctot+1:ctot+c) = wi;
                stressrv(ctot+1:ctot+c,:) = stressrvi;
                stressR(ctot+1:ctot+c,:) = stressRi;
            end
        else                                        % Tall arrays  
            if stresstype == 1
                fprintf(f.stressR,[fw,' ',fw,' ',fw,' ',fw,' ',fw,' ',fw,'\n'],stressRi');
                fprintf(f.fzH,['%s','\n'],fzHi{:});
                fprintf(f.w,['%d','\n'],wi');
            elseif stresstype == 2
                fprintf(f.philoadaxisB,[fw,'\n'],phii');
                fprintf(f.thloadaxisB,[fw,'\n'],thi');
                fprintf(f.fzloadaxis,['%s','\n'],fzloadaxisi{:});
                fprintf(f.w,['%d','\n'],wi');
                fprintf(f.nloadaxis,[fw,' ',fw,' ',fw,'\n'],nloadaxisi');
                fprintf(f.stressR,[fw,' ',fw,' ',fw,' ',fw,' ',fw,' ',fw,'\n'],stressRi');                
            else
                fprintf(f.fzstress,['%s','\n'],fzstressi{:});
                fprintf(f.w,['%d','\n'],wi');
                fprintf(f.stressrv,[fw,' ',fw,' ',fw,'\n'],stressrvi');
                fprintf(f.stressR,[fw,' ',fw,' ',fw,' ',fw,' ',fw,' ',fw,'\n'],stressRi');
            end
        end
        ctot = ctot+c;                              % Counter
    end
    
    %%% Grouping all tests data
    if ~txt                                 % Normal arrays
        if stresstype == 1
            fzH(ctot+1:end) = [];
            w(ctot+1:end) = [];    
            stressR(ctot+1:end,:) = [];
        elseif stresstype == 2
            philoadaxisB(ctot+1:end) = [];
            thloadaxisB(ctot+1:end) = [];
            fzloadaxis(ctot+1:end) = [];
            w(ctot+1:end) = [];    
            nloadaxis(ctot+1:end,:) = [];
            stressR(ctot+1:end,:) = [];
        else
            fzstress(ctot+1:end) = [];
            w(ctot+1:end) = [];    
            stressrv(ctot+1:end,:) = [];
            stressR(ctot+1:end,:) = [];
        end
    else
        ffields = fieldnames(f);
        for i = 1:length(ffields)
            fclose(f.(ffields{i}));
        end
        if stresstype == 1
            fzH = readcell('temp/fzH.txt');
            w = readmatrix('temp/w.txt');
            stressR = readmatrix('temp/stressR.txt');
        elseif stresstype == 2
            philoadaxisB = readmatrix('temp/philoadaxisB.txt');
            thloadaxisB = readmatrix('temp/thloadaxisB.txt');
            fzloadaxis = readcell('temp/fzloadaxis.txt');
            w = readmatrix('temp/w.txt');
            nloadaxis = readmatrix('temp/nloadaxis.txt');
            stressR = readmatrix('temp/stressR.txt');
        else
            fzstress = readcell('temp/fzstress.txt');
            w = readmatrix('temp/w.txt');
            stressrv = readmatrix('temp/stressrv.txt');
            stressR = readmatrix('temp/stressR.txt');
        end
    end
    if stresstype == 2
        if size(philoadaxisB,1) < size(philoadaxisB,2)              % Rearrange if needed
            philoadaxisB = philoadaxisB';
            thloadaxisB = thloadaxisB';
            fzloadaxis = fzloadaxis';
            w = w';        
        end
    elseif stresstype == 3
        if size(stressrv,1) < size(stressrv,2)
            stressrv = stressrv';
            fzstress = fzstress';
            w = w';        
        end
    end


    %%% Output: Complete set of test configurations
                                                                        % IDs, fzs and weights
    tests.id = [replicatemat(gb.id,Nloads),...                          % [Disorientation id, grain boundary id, stress id]
        cell2mat(arrayfun(@(x) 1:x,Nloads','UniformOutput',false))'];
    tests.fzdichromatic = replicatecell(gb.fzdichromatic,Nloads);
    tests.fzboundary = replicatecell(gb.fzboundary,Nloads);
    if stresstype == 1
        tests.fzH = fzH;
    elseif stresstype == 2
        tests.fzloadaxis = fzloadaxis;
    elseif stresstype == 3
        tests.fzstress = fzstress;
    end
    tests.symdichromatic = replicatecell(gb.symdichromatic,Nloads);  
    tests.symgboundary = replicatecell(symgboundary,Nloads);          
    tests.w = [replicatemat(gb.w,Nloads),w];
    
    tests.euleranglesI = replicatemat(gb.euleranglesI,Nloads);          % dis fields
    tests.euleranglesII = replicatemat(gb.euleranglesII,Nloads);
    tests.l = replicatemat(gb.l,Nloads);
    tests.th = replicatemat(gb.th,Nloads);
    tests.rv = replicatemat(gb.rv,Nloads);

    tests.thrangeboundary = replicatemat(gb.thrangeboundary,Nloads);    % gb fields
    tests.symx = replicatemat(gb.symx,Nloads);
    tests.symy = replicatemat(gb.symy,Nloads);
    tests.symz = replicatemat(gb.symz,Nloads);
    tests.symx2 = replicatemat(gb.symx2,Nloads);
    tests.phiboundary = replicatemat(gb.phiboundary,Nloads);
    tests.thboundary = replicatemat(gb.thboundary,Nloads);
    tests.nboundary = replicatemat(gb.nboundary,Nloads);
    tests.GBangles = replicatemat(gb.GBangles,Nloads);
    tests.Bx = replicatemat(gb.Bx,Nloads);
    tests.By = replicatemat(gb.By,Nloads);
    tests.Bz = replicatemat(gb.Bz,Nloads);
    tests.Bx2 = replicatemat(gb.Bx2,Nloads);
    tests.thrangeloadaxis = replicatemat(thrangeloadaxis,Nloads);  
    
    if stresstype == 2                                                  % test fields
        tests.philoadaxisB = philoadaxisB;
        tests.thloadaxisB = thloadaxisB;
        tests.nloadaxis = nloadaxis;
    elseif stresstype == 3
        tests.stressrv = stressrv;
    end
    tests.stressR = stressR;
end


% Functions to replicate variables
function m_new = replicatemat(m,n)
    m_new = cell2mat(arrayfun(@(r) repmat(m(r,:),n(r),1),(1:numel(n))','UniformOutput',false));
end

function c_new = replicatecell(c,n)
    uc = unique(c);
    uids = cell2mat(cellfun(@(x) find(strcmp(x,uc)),c,'UniformOutput',false));
    udisrep = replicatemat(uids,n);
    c_new = arrayfun(@(x) {uc{x}},udisrep);
end