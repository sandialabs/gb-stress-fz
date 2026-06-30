%%% Example - (disorientation + boundary plane) grid

% Inputs
steps_dis = 3;                              % Steps in disorientation FZ
steps_bp = 8;                               % Steps in boundary plane FZ

% Grids
dis = grid_disorientation(steps_dis);       % Disorientation 
gb = grid_boundaryplane(dis,steps_bp);      % Boundary plane
% - Other options:
% gb = grid_boundaryplane(dis,'unique');
% gb = grid_boundaryplane(dis,steps_bp,'symmetry',false);   % Boundary plane - without considering symmetry

% Plots
plot_disorientation(dis,'LabelRegion',true);% Disorientation plot

plotpoints = 'A';                           % Plot for 'A' configuration
plot_boundaryplane(gb,plotpoints);
plot_boundaryplane(gb,plotpoints,'PlotType','3D');

plotpoints = 'OE';                          % Plot for all 'OE' configurations
plot_boundaryplane(gb,plotpoints);

plotpoints = {'OA',1};                      % Plot for only the first 'OA' configuration
plot_boundaryplane(gb,plotpoints);

% For easier data visualization: struct2table(gb)
