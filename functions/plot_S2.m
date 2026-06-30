%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot an S2 fundamental zone (FZ).
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Inputs:
% plottype          '2D' projection, or '3D' patch of a sphere
% th                Polar angles of each point
% phi               Azimuthal angles of each point
% results           Property to color the FZ points according to. 
%                   nan vector to ignore color
% thrange           Polar angle covered by the FZ
% xlbl,ylbl,zlbl,x22lbl     Axes labels
% ttl               Figure title
%
%%% Outputs:
% fig               Cell with the resulting figures
%

function fig = plot_S2(plottype,th,phi,results,thrange,xlbl,ylbl,zlbl,x2lbl,ttl)

    assert(strcmp(plottype,'2D') || strcmp(plottype,'3D'),"plottype must be '2D' or '3D'")
    assert(size(th,1)==size(phi,1) && size(th,1)==size(results,1) &&...
           size(th,2)==1 && size(phi,2)==1 && size(results,2)==1,'th, phi, results must all be column vectors of the same size')
    assert(ischar(xlbl) && ischar(ylbl) && ischar(zlbl) && ischar(x2lbl) && ischar(ttl),'xlbl, ylbl, zlbl, x2lbl, ttl must all be strings')
    warning('off','MATLAB:handle_graphics:exceptions:SceneNode');

    % 2D plot
    if strcmp(plottype,'2D')

        fig = figure;
        if all(isnan(results))  % Without results vector
            polarscatter(th,sin(phi),40,'filled')  
        else                    % With results vector
            polarscatter(flipud(th),...               % Flip to plot from the edges to the center of the circle
                         flipud(sin(phi)),40,...
                         flipud(results),'filled') 
        end
        ga = gca;
        if (thrange-pi)<=10*eps
            ga.ThetaLim = [0,thrange+0.000000000001]*180/pi;
        else
        end
        thticks = thetaticks;
        if ~(abs(thticks(end)-thrange*180/pi) <= 10*eps)
            thticks = 0:(thrange*180/pi/3):thrange*180/pi;
            thetaticks(thticks)
        end
        rtcks = rticks;
        if abs(thrange-2*pi)<=10*eps                % Whole circle
            rticklabels([])
            thticklbls = cell(size(thticks'));
            thticklbls{1} = ['$',xlbl,'$'];
            thticklbls{thticks == 90} = ['$',ylbl,'$'];
            thetaticklabels(thticklbls)
            ga.TickLabelInterpreter = 'latex';
        elseif abs(thrange-pi)<=10*eps              % Half circle
            rticklbls = cell(size(rtcks'));
            rticklbls{1} = ['$',zlbl,'$'];
            rticklabels(rticklbls)
            thticklbls = cell(size(thticks'));
            thticklbls{1} = ['$',xlbl,'$'];
            thticklbls{thticks == 90} = ['$',ylbl,'$'];
            thetaticklabels(thticklbls)
            ga.TickLabelInterpreter = 'latex';
        else                                        % Circle section smaller than pi
            rticklbls = cell(size(rtcks'));
            rticklbls{1} = ['$',zlbl,'$'];
            rticklabels(rticklbls)
            thticklbls = cell(size(thticks'));
            thticklbls{1} = ['$',xlbl,'$'];
            thticklbls{end} = ['$',x2lbl,'$'];
            thetaticklabels(thticklbls)
            ga.TickLabelInterpreter = 'latex';
        end
        T = title(ttl,'FontSize',14,'Interpreter','latex');
        if abs(thrange-pi/2)<=10*eps
            ga.TitleHorizontalAlignment = 'center';
        else
            ga.TitleHorizontalAlignment = 'left';
            T.HorizontalAlignment = 'left';
        end
        ga.LineWidth = 0.75;
        set(gca,'FontName','Arial','FontSize',12)
        set(gcf,'color','w');
       
    % 3D plot
    elseif strcmp(plottype,'3D')
        fig = figure;
        if all(isnan(results))  % Without results vector
            scatter3(sin(phi).*cos(th),...
                     sin(phi).*sin(th),...
                     cos(phi),50,'filled')
        else                    % With results vector
            scatter3(sin(phi).*cos(th),...
                     sin(phi).*sin(th),...
                     cos(phi),50,...
                     results,'filled')
        end
        axis equal
        if thrange >= pi-10*eps
            xticks(-1:1:1)
            yticks(-1:1:1)
            zticks(-1:1:1)
        else
            xticks(-1:0.5:1)
            yticks(0:max(sin(th))/2:max(sin(th)))
            zticks(-1:0.5:1)
        end
    %     xlabel(['$\mathsf{x: ',xlbl,'}$'],'Interpreter','latex')    % For older versions
    %     ylabel(['$\mathsf{y: ',ylbl,'}$'],'Interpreter','latex')
    %     zlabel(['$\mathsf{z: ',zlbl,'}$'],'Interpreter','latex')
        xlabel(['$x: ',xlbl,'$'],'Interpreter','latex')
        ylabel(['$y: ',ylbl,'$'],'Interpreter','latex')
        zlabel(['$z: ',zlbl,'$'],'Interpreter','latex')

        title(ttl,'FontSize',14,'Interpreter','latex');
        set(gca,'FontName','Arial','FontSize',12)
        set(gcf,'color','w');
    end
    warning('off','MATLAB:handle_graphics:exceptions:SceneNode');
end

