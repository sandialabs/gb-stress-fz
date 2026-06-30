%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the region within an S2 fundamental zone a point 
% belongs to.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Inputs:
%   phi, th         Azimuthal and polar angles of the point
%   phimax, thmax   Maximum azimuthal and polar angles of the FZ
%   symparent       Point group of the fundamental zone of the higher level
%
%%% Outputs: 
%   fzboundary      Region within the FZ:
%                   - Points:                 
%                       O: // disorientation axis
%                       A: ⊥ disorientation axis, and // symx
%                       B: ⊥ disorientation axis, and // symx2
%                   - Edges: OA, OB, AB
%                   - interior
%   w               Weights (i.e. number of symmetric points within the FZ)
%

function [fzboundary,w] = regionS2(phi,th,phimax,thmax,symparent)
    tol = 20*eps;                        % Tolerance
    
    % Symmetry
    if thmax < 2*pi                 % Symmetry > Ci                  
        if abs(phi)<=tol
            fzboundary = 'O';           % Origin
        elseif abs(phi-phimax)<=tol && abs(th)<=tol
            fzboundary = 'A';           % Vertex 1
        elseif abs(phi-phimax)<=tol && abs(th-thmax)<=tol
            fzboundary = 'B';           % Vertex 2
        elseif abs(th)<=tol
            fzboundary = 'OA';          % Edge 1
        elseif abs(th-thmax)<=tol
            fzboundary = 'OB';          % Edge 2
        elseif abs(phi-phimax)<=tol
            fzboundary = 'AB';          % Outer edge 
        else
            fzboundary = 'interior';    % Interior
        end
    else                            % Ci symmetry
        fzboundary = 'nan';
        w = 2;
    end
    
    % Weights
    if strcmp(fzboundary,'O')
        w = 2;
    elseif strcmp(fzboundary,'A') || strcmp(fzboundary,'B')
        if strcmp(symparent,'D3d')
            w = 6;
        elseif strcmp(symparent,'C2h') || strcmp(symparent,'C2h*')
            w = 2;
        else
            w = round(pi/thmax);
        end
    elseif strcmp(fzboundary,'OA') || strcmp(fzboundary,'OB')
        if strcmp(symparent,'C2h')
            w = 4;
        else             
            w = 2*round(pi/thmax);
        end
    elseif strcmp(fzboundary,'AB')
        if strcmp(symparent,'D3d')
            w = 12;
        else
            w = 2*round(pi/thmax);
        end
    elseif strcmp(fzboundary,'interior')
        w = 4*round(pi/thmax);
    elseif strcmp(fzboundary,'nan') % Ci symmetry
        w = 2;
    end
end

