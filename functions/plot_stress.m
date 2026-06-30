%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot grain bounday-stress fundamental zones (FZ).
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Actual plotting performed by plot_S2.m for an axisymmetric stress, and 
% by plot_SO3.m for a triaxial stress.
%
%%% Inputs:
% tests             stressed grain boundary structure
% id                Points within gb to plot:
%                   'all' -> All points
%                   string -> points belonging to that disorientation region
%                   cell {FZdis,id} -> disorientation regions to plot
%                   cell {FZdis,FZgb,id} -> grain boundary regions to plot
% *Optional inputs:
% results           Property to color the FZ points according to
% 'plottype'        '2D' (default) or '3D'
% 'millerindices'   Approximate Miller indices, true (default) or false
% 'plotfz'          Plot FZ, true (default) or false. Only for triaxial stresses.
% 'labelregion'     Label each point by its region in the FZ, true or false
%                   (default). Only for triaxial stresses.
% 
%%% Outputs:
% fig               Cell with the resulting figures
%

function fig = plot_stress(tests,id,varargin)

    % Asserts
    narginchk(2,inf)
    assert(all(strcmp(id,'all')) || iscell(id) || ischar(id), "id must be 'all', a FZdis string, or cells {FZdis,id} or {FZdis,FZgb,id}")
    if iscell(id)
        assert(size(id,2) == 2 || size(id,2) == 3, 'if id is a cell then it must have 2 {FZdis,id} or 3 {FZdis,FZgb,id} elements')
    end

    % Check if there is a results input
    if mod(nargin,2) == 1
        results = varargin{1};
        assert(size(results,1)==size(tests.id,1) && size(results,2)==1,'results must be a column vector with tests size')
        resultsvector = 1;
    else 
        results = nan(size(tests.fzdichromatic,1),1);
        resultsvector = 0;
    end

    % Checck stress type
    if isfield(tests,'fzH')
        error('Use plot_boundaryplane() to plot hydrostatic stresses')
    elseif isfield(tests,'fzloadaxis')
        stresstype = 2;
    elseif isfield(tests,'fzstress')
        stresstype = 3;
    else
        error('tests must be an axisymmetric or triaxial stress structure')
    end

    % Additional options
    plottype = '2D';                        % Defaults
    millerindices = true;
    plotfz = true;
    labelregion = false;
    for i = (1+resultsvector):2:length(varargin)                
        propName = varargin{i};
        propValue = varargin{i+1};

        if stresstype == 2                      % Axisymmetric
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
                    warning('Unknown property for axisymmetric stresses: %s', propName);
            end

        else                                    % Triaxial stresses
            switch lower(propName)
    
                % Approximate crystal directions by Miller indices?
                case 'millerindices'            
                    if propValue == true
                        millerindices = true;
                    elseif propValue == false
                        millerindices = false;
                    else
                        error("PlotFZ must be true or false")
                    end
    
                % Plot the fundamental zone?
                case 'plotfz'                   
                    if propValue == true
                        plotfz = true;
                    elseif propValue == false
                        plotfz = false;
                    else
                        error("PlotFZ must be true or false")
                    end

                % Label each point according to its region in the FZ
                case 'labelregion'                   
                    if propValue == true
                        labelregion = tests.fzstress;
                    elseif propValue == false
                        labelregion = false;
                    else
                        error("LabelRegion must be false (default) or true ")
                    end

                otherwise
                    warning('Unknown property for triaxial stresses: %s', propName);
            end
        end

    end


    % Points
    if strcmp(id,'all')                 % Plot all (careful, many plots may be generated)
        id = unique(tests.id(:,[1,2]),'rows');

    elseif ischar(id)                   % Plot stresses for a given disorientation
        idu = unique(tests.id(strcmp(tests.fzdichromatic,id),[1,2]),'rows');
        if isempty(idu)
            disp(['plot_stress_axisymmetric: No "',id,'" points.'])
        end
        id = idu;

    elseif iscell(id)                   
        ids = unique(tests.id(strcmp(tests.fzdichromatic,id{1,1}) & strcmp(tests.fzboundary,id{1,2}),[1,2]),'rows');

        if size(id,2)==2                % Plot stresses for a given disorientation and boundary plane
            id = ids;

        elseif size(id,2)==3            % Plot one stress fundamental zone for a given disorientation and boundary plane
            if strcmp(id{3},'all')      % All such points
                id = ids;
            elseif id{3} > size(ids,1)
                disp(['plot_stress_axisymmetric: Not enough "',id{1},'-',id{2},'" points.'])
                id = [];
            else                        % Only one of those points
                id = ids(id{1,3},:);
            end
        end
    end

    %%%%%%%%%% Plots
    fig = cell(size(id,1),1);
    for i = 1:size(id,1)

        p = all(tests.id(:,[1,2])==id(i,:),2);      % Points to plot
        pu = find(p,1);                             % Information of the points to plot

        % Plot information
        if millerindices            % xyx axes
            ldis = trans.direction2Miller(tests.l(pu,:),1);
            ndis = trans.direction2Miller(tests.nboundary(pu,:),1);
            xlbl = trans.direction2Miller(tests.Bx(pu,:),1);
            ylbl = trans.direction2Miller(tests.By(pu,:),1);
            zlbl = trans.direction2Miller(tests.Bz(pu,:),1);
            x2lbl = trans.direction2Miller(tests.Bx2(pu,:),1);
        else
            ldis = ['[',regexprep(num2str(round(tests.l(pu,:),2)),' +',','),']'];
            ndis = ['[',regexprep(num2str(round(tests.nboundary(pu,:),2)),' +',','),']'];
            xlbl = ['[',regexprep(num2str(round(tests.Bx(pu,:),2)),' +',','),']'];
            ylbl = ['[',regexprep(num2str(round(tests.By(pu,:),2)),' +',','),']'];
            zlbl = ['[',regexprep(num2str(round(tests.Bz(pu,:),2)),' +',','),']'];
            x2lbl = ['[',regexprep(num2str(round(tests.Bx2(pu,:),2)),' +',','),']'];
        end
        ttl = ['$',tests.fzdichromatic{pu},': ',ldis,'\,',num2str(round(tests.th(pu)*180/pi,2)),'^{\circ};\,',...
                           tests.fzboundary{pu}   ,': ',ndis,'$']; % Title

        %%% Axisymmetric stresses
        if stresstype == 2          
    
            % Stress type specific information
            if max(tests.thrangeloadaxis(p)) > tests.thrangeloadaxis(pu)       % theta range
                thrange = 2*pi;
            else
                thrange = tests.thrangeloadaxis(pu);
            end
    
            % Plot
            fig{i} = plot_S2(plottype,tests.thloadaxisB(p),tests.philoadaxisB(p),results(p),thrange,xlbl,ylbl,zlbl,x2lbl,ttl);            
   
        %%% Triaxial stresses
        else                                  
            
            % Plot
            if islogical(labelregion)
                labelregionp = labelregion;
            else
                labelregionp = labelregion(p);
            end
            fig{i} = plot_SO3(tests.stressrv(p,:),results(p),tests.symgboundary{pu},'D2h',xlbl,ylbl,zlbl,ttl,'plotfz',plotfz,'labelregion',labelregionp);            
    
        end
    end

end

