%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to create boundary plane grids of Oh disorientations.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Inputs:
% dis   - Disorientations (structure from grid_disorientation)
%   euleranglesI    Euler angles of grain I [0,0,0]
%   euleranglesII   Euler angles of grain II (ZXZ)
%   l               Disorientation axis
%   th              Disorientation angle
%   rv              Rodrigues vector
%   fzdichromatic   Location in the disorientation fundamental zone
% steps             Number of boundary plane steps:
%                     single integer - quasi-equispaced grid 
%                     'unique' - One point for each region of the FZ
% *Optional inputs:
% 'symmetry'        Consider crystal symmetry to create the grid?
%                   false -> no; 2*pi range
%                   true (default) -> yes; only fundamental zone
%                   2 -> yes; only fundamental zone but do not delete
%                        edge points that appear multiple times
%                   
%%% Outputs:
% gb                Grain boundary structure:
%   id              Unique id for each condition [id_dis, id_gb]
%   []              All inputs from dis
%   fzboundary      Location in the boundary plane FZ
%   symgboundary    Point group of the grain boundary
%   symx,symy,symz,symx2    Axes in the symmetry coordinate system
%   thrangeboundary         Angle range of the boundary plane FZ
%   phiboundary,thboundary  Spherical coordinates (phi,th) of the boundary planes in the symmetry coordinate system
%   nboundary       Miller indices of the normal to the boundary plane in grain I
%   GBangles        Angles [thB,phiB] of the boundary normal
%   Bx,By,Bz,Bx2    Symmetry axes of the grain boundary plane coordinate system
%   thrangeloadaxis Angle range of the load axis FZ
%

function gb = grid_boundaryplane(varargin)

    [dis,steps] = varargin{1:2};
    % Defaults
    sym = 1;
    % Additional options
    if length(varargin) > 2
        for i = 3:2:length(varargin)      
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

                otherwise
                    error('processOptions:UnrecognizedOption', ...
                          'Unrecognized option: ''%s''', optionName);
            end
        end
    end

    assert(isstruct(dis) && isfield(dis,'l') && isfield(dis,'th') && isfield(dis,'fzdichromatic'),'dis must be a structure with specific fields (see grid_disorientation.m)')
    assert((~ischar(steps) && length(steps)==1 && floor(steps)==steps) || ...
           strcmp(steps,'unique'), ...
           "steps must be a single integer or 'unique'")
    
    lsize = size(dis.l,1);                                  % Number of disorientations
    if ~ischar(steps)                                       % Maximum number of steps
        nsteps = steps;
    else
        nsteps = 3;
    end

    % Boundary plane symmetry range
    if sym >= 1                                             % thrangeboundary 
        thrangeboundary = dis.thrangeboundary;              % Using disorientation symmetry
    else                                                    % Whole sircle
        thrangeboundary = zeros(size(dis.thrangeboundary)) + 2*pi;
    end

    % Boundary plane grid along the symmetry axes
    Nboundaries = nan(lsize,1);
    [phiboundary,thboundary,w] = deal(nan(lsize*4*max(nsteps,1)^2,1));
    fzboundary = cell(lsize*4*max(nsteps,1)^2,1);
    nboundary = nan(lsize*4*max(nsteps,1)^2,3);
    GBangles = nan(lsize*4*max(nsteps,1)^2,2);
    ctot = 0;   % Counter
        
    % Loop for each misorientation
    for i = 1:lsize  

        %%% Grid
        if ~ischar(steps)           % CASE: quasi-equispaced grid
            phimax = pi/2;                                  % Grid range
            thmax = thrangeboundary(i);                     % Range of the boundary normal in S^2
            [phii,thi,c] = gridS2(phimax,thmax,nsteps,[]);  % Grid over S^2                    
        else                        % CASE: unique
            phimax = pi/2;
            thmax = dis.thrangeboundary(i);     
            [phii,thi,c] = gridS2(phimax,thmax,'unique',dis.symdichromatic{i}); % Grid over S^2                    
        end
        if all(thi==0), phii = 0; thi = 0; c = 1; end
        if isnan(thi), thi = 0; end

        %%% Nature of the point in the boundary plane fundamental zone
        fzboundaryi = cell(size(phii));
        wi = nan(size(phii));
        for j = 1:length(phii)
            [fzboundaryi{j},wi(j)] = regionS2(phii(j),thi(j),phimax,thmax,dis.symdichromatic{i});
        end

        %%% Boundary plane normal
        nboundarysymi= [sin(phii).*cos(thi),sin(phii).*sin(thi),cos(phii)];   % In sym coordinate system
        nboundaryi = nan(length(phii),3);
        GBanglesi = nan(length(phii),2);

        rsym2I = dis.symx(i,:)' * [1,0,0] + ...                 % Rotation matrix (symmetry -> grain I)
                 dis.symy(i,:)' * [0,1,0] + ...
                 dis.symz(i,:)' * [0,0,1];
        for j = 1:length(phii)
            nboundaryi(j,:) = (rsym2I * nboundarysymi(j,:)')';  % Boundary plane normal in grain I coordinate system ...
            [~,thB,phiB] = trans.cart2sph(nboundaryi(j,:));  % ... and its corresponding GBangles [thB,phiB]
            GBanglesi(j,:) = [thB,phiB];
        end
        nboundary(ctot+1:ctot+c,:) = nboundaryi;            % Saving values
        GBangles(ctot+1:ctot+c,:) = GBanglesi;

        %%% Conglomerating and saving
        Nboundaries(i) = c;                                 % Number of boundary planes per disorientation
        phiboundary(ctot+1:ctot+c) = phii;                  % Saving values
        thboundary(ctot+1:ctot+c) = thi;
        fzboundary(ctot+1:ctot+c) = fzboundaryi;
        w(ctot+1:ctot+c) = wi;
        ctot = ctot+c;                                      % Counter
    end
    phiboundary(ctot+1:end) = [];
    thboundary(ctot+1:end) = [];
    fzboundary(ctot+1:end) = [];
    w(ctot+1:end) = [];
    nboundary(ctot+1:end,:) = [];
    GBangles(ctot+1:end,:) = [];
    if size(phiboundary,1) < size(phiboundary,2)
        phiboundary = phiboundary';
        thboundary = thboundary';
    end

    %%% Miller indices of the boundary normal in grain I cs
    fzdichromatic = replicatecell(dis.fzdichromatic,Nboundaries);
    symdichromatic = replicatecell(dis.symdichromatic,Nboundaries);
    thrangeboundary = replicatemat(thrangeboundary,Nboundaries);
    symx = replicatemat(dis.symx,Nboundaries);
    symy = replicatemat(dis.symy,Nboundaries);
    symz = replicatemat(dis.symz,Nboundaries);
    symx2 = replicatemat(dis.symx2,Nboundaries);

    symgboundary = cell(size(nboundary,1),1);        % Symmetry
    thrangeloadaxis = nan(size(nboundary,1),1);      % th range for the boundary FZ
    [Bx,By,Bz,Bx2] = deal(nan(size(nboundary,1),3)); % Symmetry axes
    for i = 1:size(nboundary,1)
        [symgboundary{i},thrangeloadaxis(i),Bx(i,:),By(i,:),Bz(i,:),Bx2(i,:)] = fccFZ.stressaxissymmetry(nboundary(i,:),...
                                                                    fzdichromatic{i},symdichromatic{i},fzboundary{i},thrangeboundary(i), ...
                                                                    symx(i,:),symy(i,:),symz(i,:),symx2(i,:),true);
    end 


    % Output: Complete set of grain boundary configurations
                                                                    % IDs, fzs and weights
    gb.id = [replicatemat((1:lsize)',Nboundaries),...               % [Disorientation id, grain boundary id]
        cell2mat(arrayfun(@(x) 1:x,Nboundaries','UniformOutput',false))'];
    gb.fzdichromatic = fzdichromatic;
    gb.fzboundary = fzboundary;
    gb.symdichromatic = symdichromatic;        
    gb.symgboundary = symgboundary;        
    gb.w = [replicatemat(dis.w,Nboundaries),w];

    gb.euleranglesI = replicatemat(dis.euleranglesI,Nboundaries);   % dis fields
    gb.euleranglesII = replicatemat(dis.euleranglesII,Nboundaries);
    gb.l = replicatemat(dis.l,Nboundaries);
    gb.th = replicatemat(dis.th,Nboundaries);
    gb.rv = replicatemat(dis.rv,Nboundaries);
    gb.symx = symx;
    gb.symy = symy;
    gb.symz = symz;
    gb.symx2 = symx2;
    gb.thrangeboundary = thrangeboundary; 

    gb.phiboundary = phiboundary;                                   % gb fields
    gb.thboundary = thboundary;
    gb.nboundary = nboundary;
    gb.GBangles = GBangles;
    gb.Bx = Bx;
    gb.By = By;
    gb.Bz = Bz;
    gb.Bx2 = Bx2;
    gb.thrangeloadaxis = thrangeloadaxis;

    
    % Deleting repeated points from FZs without a symmetry plane 
    % perpendicular to the main symmetry axis (z) have 
    if sym == 1
        del = zeros(size(gb.fzdichromatic));
        % C2h symmetry, values along theta=180 deg (OB & B)
        % -> equivalent point along theta=0 deg (OA & A)
        del = del | strcmp(gb.symdichromatic,'C2h') & ...
                    (strcmp(gb.fzboundary,'OB') | strcmp(gb.fzboundary,'B'));

        % D3d symmetry (OE points), values along AB for theta>30 deg
        % -> equivalent point at AB for theta < 30 deg
        del = del | strcmp(gb.symdichromatic,'D3d') & ...
                    ((strcmp(gb.fzboundary,'B')) | ...
                     (strcmp(gb.fzboundary,'AB') & thboundary>pi/6));
        
        % Deleting the points
        gbfields = fieldnames(gb);
        for i = 1:length(gbfields)
            gb.(gbfields{i})(del,:) = [];
        end
    else    % includeedges -> the weights are no longer accurate
        gb.w = nan(size(gb.w));
    end

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


