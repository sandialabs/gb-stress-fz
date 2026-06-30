%%% Example - disorientation grid

% Inputs
steps = 8;                              % Steps in disorientation FZ    

% Grid
dis = grid_disorientation(steps);      
% - Other options:
% dis = grid_disorientation('twin'); 
% dis = grid_disorientation('unique'); 

% Plots
plot_disorientation(dis);
plot_disorientation(dis,'LabelRegion',true);

% For easier data visualization: struct2table(dis)
