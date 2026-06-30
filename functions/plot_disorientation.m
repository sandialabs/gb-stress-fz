%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot disorientation fundamental zones (FZ).
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Actual plotting performed by plot_SO3.m
%
%%% Inputs:
% dis               disorientation structure
% *Optional inputs:
% results           Property to color the FZ points according to
% 'millerindices'   Approximate Miller indices, true (default) or false
% 'plotfz'          Plot FZ, true (default) or false
% 'labelregion'     Label each point by its region in the FZ, true or false
%                   (default)
% 
%%% Outputs:
% fig               Resulting figure
%

function fig = plot_disorientation(dis,varargin)
    
    % Asserts
    narginchk(1,inf)

    % Check if there is a results input
    if mod(nargin,2) == 0
        results = varargin{1};
        assert(size(results,1)==size(dis.fzdichromatic,1) && size(results,2)==1,'results must be a column vector with dis size')
        resultsvector = 1;
    else 
        results = nan(size(dis.fzdichromatic,1),1);
        resultsvector = 0;
    end

    % Additional options
    millerindices = true;                       % Defaults
    plotfz = true;
    labelregion = false;
    for i = (1+resultsvector):2:length(varargin)                
        propName = varargin{i};
        propValue = varargin{i+1};

        switch lower(propName)

            % Approximate crystal directions by Miller indices?
            case 'millerindices'            
                if propValue == true
                    millerindices = true;
                elseif propValue == false
                    millerindices = false;
                else
                    error("MillerIndices must be true (default) or false")
                end

            % Plot the fundamental zone?
            case 'plotfz'                   
                if propValue == true
                    plotfz = true;
                elseif propValue == false
                    plotfz = false;
                else
                    error("PlotFZ must be true (default) or false")
                end

            % Label each point according to its region in the FZ
            case 'labelregion'                   
                if propValue == true
                    labelregion = dis.fzdichromatic;
                elseif propValue == false
                    labelregion = false;
                else
                    error("LabelRegion must be false (default) or true ")
                end

            otherwise
                warning('Unknown property for triaxial stresses: %s', propName);
        end
    end

    %%%%%%%%%% Plots

    % Plot information
    if millerindices            % xyx axes
        xlbl = trans.direction2Miller([1,0,0],1);
        ylbl = trans.direction2Miller([0,1,0],1);
        zlbl = trans.direction2Miller([0,0,1],1);
    else
        xlbl = ['[',regexprep(num2str([1,0,0]),' +',','),']'];
        ylbl = ['[',regexprep(num2str([0,1,0]),' +',','),']'];
        zlbl = ['[',regexprep(num2str([0,0,1]),' +',','),']'];
    end
    ttl = ''; % Title

    % Plot
    fig = plot_SO3(dis.rv,results,'Oh','Oh',xlbl,ylbl,zlbl,ttl,'plotfz',plotfz,'labelregion',labelregion); 

end