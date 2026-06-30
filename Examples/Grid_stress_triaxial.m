%%% Example - (disorientation + boundary plane + triaxial stress) grid

% Inputs                         
type_load = 'triaxial';                     % Triaxial stress (Try also 'shear' or 'tensor')
stress = [3,1,2]';                          % Magnitude of the eigenstresses
steps_dis = 3;                              % Steps in disorientation FZ
steps_bp = 4;                               % Steps in boundary plane FZ
steps_load = 6;                             % Steps in stress FZ

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
plot_stress(tests,plotpoints,'LabelRegion',true);

plotpoints = 'E';                                       % Plot for all disorientation 'E' configurations
plot_stress(tests,plotpoints);

% For easier data visualization: struct2table(tests)