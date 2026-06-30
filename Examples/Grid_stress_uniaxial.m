%%% Example - (disorientation + boundary plane + axisymmetric stress) grid

% Inputs                         
type_load = 'uniaxial';                     % Uniaxial stress (Try also 'biaxial')
stress = 1;                                 % Magnitude of the uniaxial stress
steps_dis = 4;                              % Steps in disorientation FZ
steps_bp = 4;                               % Steps in boundary plane FZ
steps_load = 8;                             % Steps in load axis FZ

% Grids
dis = grid_disorientation(steps_dis);                   % Disorientation 
gb = grid_boundaryplane(dis,steps_bp);                  % Boundary plane
tests = grid_stress(gb,type_load,stress,steps_load);    % Stress
% - Other options:
% tests = grid_stress(gb,type_load,stress,steps_load,'symmetry',false);   % Stress - without considering symmetry

% Plots
plot_disorientation(dis,'LabelRegion',true);            % Disorientation plot
plot_boundaryplane(gb,'E');                             % Boundary plane plot

plotpoints = {'A','OA',1};                              % Plot for the first {disorientation 'A', boundary plane 'OA'} configurations
plot_stress(tests,plotpoints);
plot_stress(tests,plotpoints,'PlotType','3D');

plotpoints = 'E';                                       % Plot for all disorientation 'E' configurations
plot_stress(tests,plotpoints);

% For easier data visualization: struct2table(tests)
