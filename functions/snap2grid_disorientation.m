%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to snap disorientations to points in a predefined grid
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Inputs:
% dis             Disorientation conditions - toFZ_disorientation()
% dis_grid        Predefined grid - grid_disorientation()
% *NOTE - Both structures must be mapped to the FZ
%
%%% Outputs:
% idgrid          ID of the closest grid point
% thgrid          Angle [deg] to that closest grid point
% dis_snapped     dis structure for each of the closest gb_grid points
%

function [idgrid,thgrid,dis_snapped] = snap2grid_disorientation(dis,dis_grid)
    %%%%% Asserts & initialization   
    assert(isstruct(dis),'dis must be a disorientation structure')
    assert(isstruct(dis_grid),'dis_grid must be a disorientation structure')

    % Number of conditions
    n_in = size(dis.euleranglesI,1);                            % Conditions
    n_grid = size(dis_grid.euleranglesI,1);                     % Grid points
    
    %%%%% Building disorientations for each point
    rdis = cell(n_in,1);
    rdisgrid = cell(n_grid,1);
    for i = 1:n_in                                              % Conditions
        rdis{i} = trans.raxisangle(dis.l(i,:),dis.th(i));
    end
    for i = 1:n_grid                                            % Grid points
        rdisgrid{i} = trans.raxisangle(dis_grid.l(i,:),dis_grid.th(i));
    end
   
    %%%%% Finding the closest grid point
    [idgrid,thgrid] = deal(nan(n_in,1));
    for i = 1:n_in
        min_th = 1000;                                          % Preallocating idgrid
        min_idgrid = 0;                
        for j = 1:n_grid                                        % Searching each grid point
            rdis2disgrid = rdis{i}*rdisgrid{j}';                % Rotation between disorientations
            [~,th] = trans.r2axisangle(rdis2disgrid);           % Angle
            if th < min_th                                      % Save minimum
                min_idgrid = j;
                min_th = th;
            end
        end
        idgrid(i) = min_idgrid;
        thgrid(i) = min_th * 180/pi;
    end

    %%%%% Remap to the disorientation FZ after snapping to grid.
    % The fundamental zone of the corresponding point in dis_grid is used, 
    % which always has the same or higher symmetry.
    dis_snapped = dis;
    disfields = fieldnames(dis_grid);
    for i = 1:length(disfields)
        dis_snapped.(disfields{i}) = dis_grid.(disfields{i})(idgrid,:);
    end
    dis_snapped.w = nan(n_in,1);
end
