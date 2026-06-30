%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to create a grid in S2 orientation space.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Inputs:
% phimax        Maximum aximuthal angle [rad]
% thmax         Maximum polar angle [rad]
% steps         Grid steps in each dimension (every 90deg rotation)
% pg            Point group symmetry  (only needed for steps='unique')     
%
%%% Outputs:
% phi           Azimuthal angles [rad]
% th            Polar angle [rad]
% c             Number of points in the grid
%

function [phi,th,c] = gridS2(phimax,thmax,steps,pg)
    assert(phimax<=pi,'phimax must be <= pi')
    assert(thmax<=2*pi,'thmax must be <= 2*pi')
    assert((~ischar(steps) && length(steps)==1 && all(floor(steps)==steps)) || ...
       strcmp(steps,'unique'), ...
       "steps must be an integer or the string 'unique'")
    
    %%%%%%%%%%%%%%%%% GRID
    if isnumeric(steps)
        % Factors in thvmax steps to preserve homogeneous sampling density
        thvfactor = thmax/(pi/2);
        if thmax < 2*pi
            thvmax = thmax;
        else
            thvmax = 2*pi-thmax/(thvfactor*steps);
        end
    
        % Sampling
        [phi,th] = deal(nan(max(1,4*steps^2),1));               % Preallocating
        stepsphi = ceil(steps*phimax/(pi/2));                   % Number of steps in phi
        phiv = 0:(phimax/stepsphi):phimax;                      % Linear sampling on phi (irregular on r2D)
        r2Dv = sin(phiv);                                       % Corresponding radius in the 2D projection
        c = 0;  % Counter
    
        for j = 1:length(phiv)
            steps_thvj = ceil(r2Dv(j)*steps);                   % Number of steps proportial to the circumference
            thv = 0:thmax/ceil(thvfactor*steps_thvj):thvmax;
            if length(thv) <= 1                                 % At least point
                thv = 0;
            end
            for k = 1:length(thv)
                c = c+1;
                phi(c,1) = phiv(j);
                th(c,1) = thv(k);                               % Symmetry range
            end
        end
        phi(c+1:end) = [];                                      % Trimming
        th(c+1:end)= [];

    %%%%%%%%%%%%%%%%% UNIQUE
    elseif(strcmp(steps,'unique'))
        if strcmp(pg,'Ci')
            c = 1;                      % No symmetry
            phi = 0;                    % interior
            th = 0;
        elseif strcmp(pg,'C2h')
            c = 3;                      % 2-fold rotation only
            phi = [0, pi/4,pi/2]';      % O, AB, interior
            th  = [0,    0,   0]';
        else
            c = 7;                      % All other symmetries
            phi = [0, pi/2,   pi/2, pi/4,  pi/4,    pi/2,    pi/4]';
            th  = [0,    0,  thmax,    0, thmax, thmax/2, thmax/2]';
        end

    end

end