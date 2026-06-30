%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symmetries of the different point group symmetries.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Proper symmetry rotations were incorporated from MTEX with a script 
% provided at the end of this file. Another script is provided to plot the
% symmetries of a given point group.
%
% Special cases:
% - Point group D8h is not crystallographic so it was added manually.
% - The configuration of point group D3d is rotated by 30deg from the 
%   standard representation.
% - Point group C2h* corresponds to C2h observed such that the mirror plane
%   appears as a horizontal line.
%
%%% Inputs
% p         Point group
%
%%% Outputs:
% r         List of symmetry rotation matrices
% n         Number of symmetry rotations
%

function [r,n] = symmetryrotations(p)

    % Crystallographic - from MTEX
    if strcmp(p,'Oh')
        n = 24;
        r{1} = [1 0 0; 0 1 0; 0 0 1];
        r{2} = [0 1 0; 0 0 1; 1 0 0];
        r{3} = [0 0 1; 1 0 0; 0 1 0];
        r{4} = [0 1 0; 1 0 0; 0 0 -1];
        r{5} = [0 0 1; 0 1 0; -1 0 0];
        r{6} = [1 0 0; 0 0 1; 0 -1 0];
        r{7} = [0 1 0; -1 0 0; 0 0 1];
        r{8} = [0 0 1; 0 -1 0; 1 0 0];
        r{9} = [1 0 0; 0 0 -1; 0 1 0];
        r{10} = [1 0 0; 0 -1 0; 0 0 -1];
        r{11} = [0 1 0; 0 0 -1; -1 0 0];
        r{12} = [0 0 1; -1 0 0; 0 -1 0];
        r{13} = [-1 0 0; 0 -1 0; 0 0 1];
        r{14} = [0 -1 0; 0 0 -1; 1 0 0];
        r{15} = [0 0 -1; -1 0 0; 0 1 0];
        r{16} = [0 -1 0; -1 0 0; 0 0 -1];
        r{17} = [0 0 -1; 0 -1 0; -1 0 0];
        r{18} = [-1 0 0; 0 0 -1; 0 -1 0];
        r{19} = [0 -1 0; 1 0 0; 0 0 1];
        r{20} = [0 0 -1; 0 1 0; 1 0 0];
        r{21} = [-1 0 0; 0 0 1; 0 1 0];
        r{22} = [-1 0 0; 0 1 0; 0 0 -1];
        r{23} = [0 -1 0; 0 0 1; -1 0 0];
        r{24} = [0 0 -1; 1 0 0; 0 -1 0];

    elseif strcmp(p,'D8h')         % Non-crystallographic - not from MTEX
        n = 16;
        c = cosd(45);
        r{1} = [1 0 0; 0 1 0; 0 0 1];
        r{2} = [0 -1 0; 1 0 0; 0 0 1];
        r{3} = [-1 0 0; 0 -1 0; 0 0 1];
        r{4} = [0 1 0; -1 0 0; 0 0 1];
        r{5} = [c -c 0; c c 0; 0 0 1];
        r{6} = [-c -c 0; c -c 0; 0 0 1];
        r{7} = [-c c 0; -c -c 0; 0 0 1];
        r{8} = [c c 0; -c c 0; 0 0 1];
        r{9} = [1 0 0; 0 -1 0; 0 0 -1];
        r{10} = [-1 0 0; 0 1 0; 0 0 -1];
        r{11} = [0 1 0; 1 0 0; 0 0 -1];
        r{12} = [0 -1 0; -1 0 0; 0 0 -1];
        r{13} = [c c 0; c -c 0; 0 0 -1];
        r{14} = [-c c 0; c c 0; 0 0 -1];
        r{15} = [-c -c 0; -c c 0; 0 0 -1];
        r{16} = [c -c 0; -c -c 0; 0 0 -1];

    elseif strcmp(p,'D6h')
        n = 12;
        c = cosd(30);
        r{1} = [1 0 0; 0 1 0; 0 0 1];
        r{2} = [0.5 -c 0; -c -0.5 0; 0 0 -1];
        r{3} = [0.5 c 0; -c 0.5 0; 0 0 1];
        r{4} = [-0.5 -c 0; -c 0.5 0; 0 0 -1];
        r{5} = [-0.5 c 0; -c -0.5 0; 0 0 1];
        r{6} = [-1 0 0; 0 1 0; 0 0 -1];
        r{7} = [-1 0 0; 0 -1 0; 0 0 1];
        r{8} = [-0.5 c 0; c 0.5 0; 0 0 -1];
        r{9} = [-0.5 -c 0; c -0.5 0; 0 0 1];
        r{10} = [0.5 c 0; c -0.5 0; 0 0 -1];
        r{11} = [0.5 -c 0; c 0.5 0; 0 0 1];
        r{12} = [1 0 0; 0 -1 0; 0 0 -1];

    elseif strcmp(p,'D4h')
        n = 8;
        r{1} = [1 0 0; 0 1 0; 0 0 1];
        r{2} = [1 0 0; 0 -1 0; 0 0 -1];
        r{3} = [0 1 0; -1 0 0; 0 0 1];
        r{4} = [0 -1 0; -1 0 0; 0 0 -1];
        r{5} = [-1 0 0; 0 -1 0; 0 0 1];
        r{6} = [-1 0 0; 0 1 0; 0 0 -1];
        r{7} = [0 -1 0; 1 0 0; 0 0 1];
        r{8} = [0 1 0; 1 0 0; 0 0 -1];

    elseif strcmp(p,'D2h')
        n = 4;
        r{1} = [1 0 0; 0 1 0; 0 0 1];
        r{2} = [1 0 0; 0 -1 0; 0 0 -1];
        r{3} = [-1 0 0; 0 -1 0; 0 0 1];
        r{4} = [-1 0 0; 0 1 0; 0 0 -1];

    elseif strcmp(p,'D3d')  % The D3d configurations are rotated from the 
                            % standard representation
        n = 6;
        % Standard representation
        % r{1} = [1 0 0; 0 1 0; 0 0 1];
        % r{2} = [-0.5 -0.86603 0; -0.86603 0.5 0; 0 0 -1];
        % r{3} = [-0.5 0.86603 0; -0.86603 -0.5 0; 0 0 1];
        % r{4} = [-0.5 0.86603 0; 0.86603 0.5 0; 0 0 -1];
        % r{5} = [-0.5 -0.86603 0; 0.86603 -0.5 0; 0 0 1];
        % r{6} = [1 0 0; 0 -1 0; 0 0 -1];
        % Rotated by 30 deg
        c = cosd(30);
        r{1} = [1 0 0; 0 1 0; 0 0 1];
        r{2} = [0.5 c 0; c -0.5 0; 0 0 -1];
        r{3} = [-0.5 c 0; -c -0.5 0; 0 0 1];
        r{4} = [0.5 -c 0; -c -0.5 0; 0 0 -1];
        r{5} = [-0.5 -c 0; c -0.5 0; 0 0 1];
        r{6} = [-1 0 0; 0 1 0; 0 0 -1];
        
    elseif strcmp(p,'C2h')
        n = 2;
        r{1} = [1 0 0; 0 1 0; 0 0 1];
        r{2} = [-1 0 0; 0 -1 0; 0 0 1];

    elseif strcmp(p,'C2h*')
        n = 2;
        r{1} = [1 0 0; 0 1 0; 0 0 1];
        r{2} = [-1 0 0; 0 1 0; 0 0 -1];

    elseif strcmp(p,'Ci')
        n = 1;
        r{1} = [1 0 0; 0 1 0; 0 0 1];


    end
end

% -------------------------------------------------------------------------
% Script to generate the proper rotations of a crystallographic point 
% group using MTEX (v. 6.0.0) 
% -------------------------------------------------------------------------
% p = 'Oh';                                       % Point group
% cs = crystalSymmetry(p);                        % MTEX symmetry object
% prot = cs.rot(find(cs.rot.i == 0));             % Only proper rotations
% n = length(prot(:));                            % Number of rotations
% disp(['        n = ',num2str(n),';'])           % PRINT n
% for i = 1:n
%     r = trans.rEuler(prot(i).phi1,...           % Rotation matrix
%                      prot(i).Phi,...
%                      prot(i).phi2); 
%     r(abs(r)<6*eps) = 0;                        % Numerical cleanup
%     rprint = ['[',num2str(r(1,:)),'; ',...      % r string
%                   num2str(r(2,:)),'; ',...
%                   num2str(r(3,:)),']'];
%     rprint = regexprep(rprint, '\s+', ' ');
%     disp(['        r{',num2str(i),'} = ',rprint,';'])   % PRINT r
% end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Script to plot the symmetry operations
% -------------------------------------------------------------------------
% v = trans.sph2cart(1,5*pi/180,80*pi/180);     % Vector
% [syms,nsym] = symmetryrotations('D4h');       % <--- select the point group
% [vsp,vsi] = deal(nan(3,nsym));
% for i = 1:nsym
%     vsp(:,i) = syms{i}*v;
%     vsi(:,i) = -syms{i}*v;
% end
% vsp = vsp';
% vsi = vsi';
% figure
% scatter3(vsp(:,1),vsp(:,2),vsp(:,3),'b','filled')
% hold on
% scatter3(vsi(:,1),vsi(:,2),vsi(:,3),'r','filled')
% axis equal
% -------------------------------------------------------------------------

