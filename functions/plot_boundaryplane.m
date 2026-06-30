%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot boundary plane fundamental zones (FZ).
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Actual plotting performed by plot_S2.m
%
%%% Inputs:
% gb                grain boundary structure
% id                Points within gb to plot:
%                   'all' -> All points
%                   string -> points belonging to that disorientation region
%                   cell {FZdis,id} -> disorientation regions to plot
% *Optional inputs:
% results           Property to color the FZ points according to
% 'plottype'        '2D' (default) or '3D'
% 'millerindices'   Approximate Miller indices, true (default) or false
% 
%%% Outputs:
% fig               Cell with the resulting figures
%

function fig = plot_boundaryplane(gb,id,varargin)
    
    % Asserts
    narginchk(2,inf)
    assert(all(strcmp(id,'all')) || iscell(id) || ischar(id), "id must be 'all', a FZdis string or cell {FZdis,id}")
    if iscell(id)
        assert(size(id,2) == 2, 'if id is a cell then it must have 2 elements: {FZdis,id}')
    end

    % Check if there is a results input
    if mod(nargin,2) == 1
        results = varargin{1};
        assert(size(results,1)==size(gb.id,1) && size(results,2)==1,'results must be a column vector with gb size')
        resultsvector = 1;
    else 
        results = nan(size(gb.fzdichromatic,1),1);
        resultsvector = 0;
    end

    % Additional options
    plottype = '2D';                        % Defaults
    millerindices = true;
    for i = (1+resultsvector):2:length(varargin)                
        propName = varargin{i};
        propValue = varargin{i+1};
        switch lower(propName)
            % 2D or 3D plot
            case 'plottype'                 
                if strcmp(propValue,'2D')
                    plottype = '2D';
                elseif strcmp(propValue,'3D')
                    plottype = '3D';
                else
                    error("PlotType must be '2D' (default) or 3D")
                end

            % Approximate crystal directions by Miller indices?
            case 'millerindices'            
                if propValue == true
                    millerindices = true;
                elseif propValue == false
                    millerindices = false;
                else
                    error("PlotFZ must be true or false")
                end

            otherwise
                warning('Unknown property: %s', propName);
        end
    end

    % Points
    if strcmp(id,'all')                     % Plot all (careful, many plots may be generated)
        id = unique(gb.id(:,1));

    elseif ischar(id)                       % Plot boundary planes for a given disorientation type
        idu = unique(gb.id(strcmp(gb.fzdichromatic,id),1));
        if isempty(idu)
            disp(['Boundary plane: No "',id,'" points.'])
        end
        id = idu;

    elseif iscell(id)                       % Plot boundary planes ofr one disorientation
        ids = unique(gb.id(strcmp(gb.fzdichromatic,id{1}),1));
        if id{2} > length(ids)
            disp(['Boundary plane: Not enough "',id{1},'" points.'])
            id = [];
        else
            id = ids(id{2});
        end
    end

    %%%%%%%%%% Plots
    fig = cell(size(id,1),1);
    for i = 1:size(id,1)
        p = gb.id(:,1)==id(i);      % Points to plot
        pu = find(p,1);             % Information of the points to plot

        % Plot information
        if millerindices            % xyx axes
            ldis = trans.direction2Miller(gb.l(pu,:),1);
            xlbl = trans.direction2Miller(gb.symx(pu,:),1);
            ylbl = trans.direction2Miller(gb.symy(pu,:),1);
            zlbl = trans.direction2Miller(gb.symz(pu,:),1);
            x2lbl = trans.direction2Miller(gb.symx2(pu,:),1);
        else
            ldis = ['[',regexprep(num2str(round(gb.l(pu,:),2)),' +',','),']'];
            xlbl = ['[',regexprep(num2str(round(gb.symx(pu,:),2)),' +',','),']'];
            ylbl = ['[',regexprep(num2str(round(gb.symy(pu,:),2)),' +',','),']'];
            zlbl = ['[',regexprep(num2str(round(gb.symz(pu,:),2)),' +',','),']'];
            x2lbl = ['[',regexprep(num2str(round(gb.symx2(pu,:),2)),' +',','),']'];
        end
        if max(gb.thboundary(p)) > gb.thrangeboundary(pu)       % theta range
            thrange = 2*pi;
        else
            thrange = gb.thrangeboundary(pu);
        end
        ttl = ['$',gb.fzdichromatic{pu},': ',ldis,'\,',num2str(round(gb.th(pu)*180/pi,2)),'^{\circ}$']; % Title

        % Plot
        fig{i} = plot_S2(plottype,gb.thboundary(p),gb.phiboundary(p),results(p),thrange,xlbl,ylbl,zlbl,x2lbl,ttl);

    end

end
