%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot an SO(3) fundamental zone (FZ) parametrized by Rodrigues
% vectors.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Inputs:
% rv                Rodrigues vector of each point
% results           Property to color the FZ points according to. 
%                   nan vector to ignore color
% pg1,pg2           Point groups of the two elements 
% xlbl,ylbl,zlbl,x22lbl     Axes labels
% ttl               Figure title
% *Optional inputs:
% 'plotfz'          Plot FZ, true (default) or false
% 'labelregion'     Label each point by its region in the FZ, true or false
%                   (default)
% 
%%% Outputs:
% fig               Cell with the resulting figures
%

function fig = plot_SO3(rv,results,pg1,pg2,xlbl,ylbl,zlbl,ttl,varargin)

    assert(size(rv,2)==3,'rv must be an nx3 array')
    assert(size(rv,1)==size(results,1),'results and rv must have the same number of rows')
    assert(ischar(xlbl) && ischar(ylbl) && ischar(zlbl) && ischar(ttl),'xlbl, ylbl, zlbl, x2lbl, ttl must all be strings')
    warning('off','MATLAB:handle_graphics:exceptions:SceneNode');


    % Additional options
    plotfz = true;
    labelregion = false;
    for i = 1:2:length(varargin)                
        propName = varargin{i};
        propValue = varargin{i+1};

        switch lower(propName)

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
                if islogical(propValue)
                    if propValue == false
                        labelregion = false;
                    end
                elseif size(propValue,1) == size(rv,1)
                    labelregion = propValue;
                else
                    error("LabelRegion must be false (default) or a ")
                end

            otherwise
                warning('Unknown property for SO(3) plot: %s', propName);
        end
    end

    % 3D plot
    [V,F] = fzSO3(pg1,pg2);         % Vertices and faces      
    fig = figure;
    if plotfz                       % Plot FZ
        for f = 1:length(F)                     
            plygn = fill3(V(F{f},1),V(F{f},2),V(F{f},3),'k','EdgeColor','k','LineWidth',1.5,'HandleVisibility','off');
            hold on                
            plygn.FaceAlpha = 0;
        end
    end
    if islogical(labelregion)       % Plot all points at once 
        if all(isnan(results))          % Without results vector
            scatter3(rv(:,1),rv(:,2),rv(:,3),20,'filled')
        else                            % With results vector
            scatter3(rv(:,1),rv(:,2),rv(:,3),20,results,'filled') 
        end
    else                            % Plot by region within the FZ
        fzu = unique(labelregion);
        [~,stringLength] = sort(cellfun(@length,fzu),'ascend');
        fzu = fzu(stringLength);
        for i = 1:length(fzu)
            clrs = distinguishable_colors(length(fzu));
            rvi = rv(cellfun(@(x) strcmp(x,fzu{i}),labelregion),:);
            scatter3(rvi(:,1),rvi(:,2),rvi(:,3),20,'o','MarkerEdgeColor',clrs(i,:),'MarkerFaceColor',clrs(i,:),'DisplayName',fzu{i})
            hold on
            legend('Location','EastOutside','FontSize',min(max(5,150/length(fzu)),12))
        end
    end
    hold on
    axis equal
%     xlabel(['$\mathsf{r_x: ',xlbl,'}$'],'Interpreter','latex')    % For older versions
%     ylabel(['$\mathsf{r_y: ',ylbl,'}$'],'Interpreter','latex')
%     zlabel(['$\mathsf{r_z: ',zlbl,'}$'],'Interpreter','latex')
    xlabel(['$r_x: ',xlbl,'$'],'Interpreter','latex')
    ylabel(['$r_y: ',ylbl,'$'],'Interpreter','latex')
    zlabel(['$r_z: ',zlbl,'$'],'Interpreter','latex')
    title(ttl,'FontSize',14,'Interpreter','latex');
    set(gca,'FontName','Arial','FontSize',12)
    set(gcf,'color','w');
    warning('off','MATLAB:handle_graphics:exceptions:SceneNode');
end