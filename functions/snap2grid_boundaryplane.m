%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to snap boundary planes to points in a predefined grid
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% We employ the geodesic octonion metric introduced by Francis et al. 
% (Francis et al, A geodesic octonion metric for grain boundaries, Acta 
% Mater. 166 (2019) 134-147). The inputs are all referenced to grain I
% and the octonion lies on the B coordinate system (z=0 is the grain 
% boundary plane).
%
%%% Inputs:
% gb              Boundary plane conditions - toFZ_boundaryplane()
% gb_grid         Predefined grid - grid_boundaryplane()
% *NOTE - gb and gb_grid structures must be mapped to the FZ
%
%%% Outputs:
% idgrid          ID of the closest grid point
% thgrid          Angle [deg] to that closest grid point
% gb_snapped      gb structure for each of the closest gb_grid points
%

function [idgrid,thgrid,gb_snapped] = snap2grid_boundaryplane(gb,gb_grid)

    % Number of conditions
    n_in = size(gb.euleranglesI,1);                             % Conditions
    n_grid = size(gb_grid.euleranglesI,1);                      % Grid points
    
    %%%%% Building grain boundary octonions for each point
    oct_in =  nan(n_in,8);
    oct_grid = nan(n_grid,8);

    % ... for each condition
    for i = 1:n_in
        % Boundary plane in the I CS
        nboundaryT = gb.nboundary(i,:)';                    % Normal to the GB
        if nboundaryT(3) < 1                                % Two vectors perpendicular to nboundary
            nboundaryT_perp1 = [ nboundaryT(2)/sqrt(nboundaryT(1)^2+nboundaryT(2)^2),...
                                -nboundaryT(1)/sqrt(nboundaryT(1)^2+nboundaryT(2)^2),...
                                 0]';                       %   _perp1 chosen perpendicular to z
        else
            nboundaryT_perp1 = [1,0,0]';
        end
        nboundaryT_perp2 = cross(nboundaryT,nboundaryT_perp1);

        % Rotation matrices: changing vectors from X to Y such that RX2Y * _X = _Y
        rI2B = [1;0;0]*nboundaryT_perp1' + ...   
               [0;1;0]*nboundaryT_perp2' + ...
               [0;0;1]*nboundaryT';

        % Grain boundary octonion
        rI  = eye(3);                                   % Disorientation in I cs
        rII = trans.raxisangle(gb.l(i,:),gb.th(i,:));
        rIB  = rI2B * rI;                               % Disorientation in B cs
        rIIB = rI2B * rII;
        qIB  = trans.r2quat(rIB);                       % Quaternions      
        qIIB = trans.r2quat(rIIB);
        oct_in(i,:) = qo.quat2oct(qIB,qIIB);            % Octonion
    end

    % ... for each grid point (same routine)
    for i = 1:n_grid
        % Boundary plane in the I CS
        nboundaryT = gb_grid.nboundary(i,:)';               % Normal to the GB
        if nboundaryT(3) < 1                                % Two vectors perpendicular to nboundary
            nboundaryT_perp1 = [ nboundaryT(2)/sqrt(nboundaryT(1)^2+nboundaryT(2)^2),...
                                -nboundaryT(1)/sqrt(nboundaryT(1)^2+nboundaryT(2)^2),...
                                 0]';                        %   _perp1 chosen perpendicular to z
        else
            nboundaryT_perp1 = [1,0,0]';
        end
        nboundaryT_perp2 = cross(nboundaryT,nboundaryT_perp1);
    
        % Rotation matrices: changing vectors from X to Y such that RX2Y * _X = _Y
        rI2B = [1;0;0]*nboundaryT_perp1' + ...   
               [0;1;0]*nboundaryT_perp2' + ...
               [0;0;1]*nboundaryT';

        % Grain boundary octonion
        rI  = eye(3);                                   % Disorientation in I cs
        rII = trans.raxisangle(gb_grid.l(i,:),gb_grid.th(i,:));
        rIB  = rI2B * rI;                               % Disorientation in B cs
        rIIB = rI2B * rII;           
        qIB  = trans.r2quat(rIB);                       % Quaternions
        qIIB = trans.r2quat(rIIB);
        oct_grid(i,:) = qo.quat2oct(qIB,qIIB);          % Octonion
    end
    
       
    %%%%% Finding the closest grid point
    % The two octonions are compared, but the grid octonion is rotated
    % along z in the range [0,2pi] rad. The minimum of that function is the
    % angle between the two grain boundary octonions.
    [idgrid,thgrid] = deal(nan(n_in,1));
    options = optimset('TolX',1e-2);                    % Tolerance (default is 1e-4)
    for i = 1:n_in
        min_th = 1000;                                          % Preallocating idgrid
        min_idgrid = 0;   
        for j = 1:n_grid                                        % Searching each grid point
            % Function for the angle rotating the grid octonion
            anglefun = @(theta) qo.Oangle(oct_in(i,:),qo.Orotz(oct_grid(j,:),theta));

            % Finding the angle
            [~,th] = fminbnd(anglefun,0,2*pi,options);          % Angle

            % Saving the minimum angle
            if th < min_th                                      % Save minimum
                min_idgrid = j;
                min_th = th;
            end
        end
        idgrid(i) = min_idgrid;
        thgrid(i) = min_th * 180/pi;
    end

    %%%%% Remap to the boundary plane FZ after snapping to grid.
    % The fundamental zone of the corresponding point in gb_grid is used, 
    % which always has the same or higher symmetry.
    gb_snapped = gb;
    gbfields = fieldnames(gb_grid);
    for i = 1:length(gbfields)
        gb_snapped.(gbfields{i}) = gb_grid.(gbfields{i})(idgrid,:);
    end
    gb_snapped.w = nan(n_in,2);
end
