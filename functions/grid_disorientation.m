%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to create disorientation grids of Oh materials.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Inputs:
% steps             Number of disorientation steps:
%                     1x2 integer vector - quasi-equispaced grid [steps_gIIrotationaxis,steps_gIIrotationangle] 
%                     'twin' - Only twinning configuration
%                     'unique' - One point for each region of the FZ
%                   
%%% Outputs:
% dis               Disorientation structure:
%   fzdichromatic   Location in the disorientation fundamental zone
%   symdichromatic  Point group of the bicrystal
%   w               Weights (times each point repeats over the whole orientation space)
%   euleranglesI    Euler angles of grain I [0,0,0]
%   euleranglesII   Euler angles of grain II (ZXZ)
%   l               Disorientation axis
%   th              Disorientation angle [rad]
%   rv              Rodrigues vector
%   symx,symy,symz,symx2    Axes in the symmetry coordinate system
%   thrangeboundary Angle range of the boundary plane FZ
%

function dis = grid_disorientation(steps)

    assert((~ischar(steps) && length(steps)==1 && all(floor(steps)==steps)) || ...
           strcmp(steps,'twin') || strcmp(steps,'unique'), ...
           "steps must be an integer, 'twin' or 'unique'")
    
    % Grain I - not rotated
    euleranglesI = [0,0,0];

    % Disorientation grid
    ge = 1;                             % Grain exchange symmetry ON
    reg = 1;                            % Output region ON
    if ~ischar(steps)                   % CASE: quasi-equispaced grid
        [rv,fzdichromatic,w] = gridSO3('Oh','Oh',steps,ge,reg);
    elseif strcmp(steps,'unique')       % CASE: unique
        [rv,fzdichromatic,w] = gridSO3('Oh','Oh',steps,ge,reg);
    elseif strcmp(steps,'twin')         % CASE: twin
        rv = [1/3,1/3,1/3];
        [fzdichromatic,w] = regionSO3('Oh','Oh',rv,ge);
    end

    % Angle-axis vectors
    l = nan(size(rv));
    th = nan(size(rv,1),1);
    for i = 1:size(rv,1)
        [l(i,:),th(i)] = trans.rodrigues2axisangle(rv(i,:));
    end
    
    % Euler angles
    euleranglesII = nan(size(l,1),3);
    for i = 1:size(l,1)
        rot = trans.raxisangle(l(i,:),th(i));                                               % Rotation matrix
        [euleranglesII(i,1),euleranglesII(i,2),euleranglesII(i,3)] = trans.r2euler(rot);    % Euler angles
    end
    
    % Symmetry and symmetry axes for each disorientation - Table 21 of [Patala 2013]
    symdichromatic = cell(size(rv,1),1);                    % Symmetry
    thrangeboundary = nan(size(rv,1),1);                    % th range for the boundary FZ
    [symx,symy,symz,symx2] = deal(nan(size(rv,1),3));       % Symmetry axes
    for i = 1:size(rv,1)
        quat = trans.rodrigues2quaternion(rv(i,:));         % Quaternion
        [symdichromatic{i},thrangeboundary(i),symx(i,:),symy(i,:),symz(i,:),symx2(i,:)] = fccFZ.boundaryplanesymmetry(quat,fzdichromatic{i},true);
    end

    % Unique disorientations
    dis.fzdichromatic = fzdichromatic;
    dis.symdichromatic = symdichromatic;
    dis.w = w;
    dis.euleranglesI = euleranglesI .* ones(size(l,1),1);
    dis.euleranglesII = euleranglesII;
    dis.l = l;
    dis.th = th;
    dis.rv = rv;
    dis.symx = symx;
    dis.symy = symy;
    dis.symz = symz;
    dis.symx2 = symx2;
    dis.thrangeboundary = thrangeboundary;

end












