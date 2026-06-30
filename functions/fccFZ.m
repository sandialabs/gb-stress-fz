%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions to assign symmetries and axes of grain boundary coordinate 
% systems.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Functions:
% boundaryplanesymmetry
% stressaxissymmetry
%

classdef fccFZ  
    methods(Static)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Symmetry axes - boundary plane
        % To determine the symmetry axes according to Table 21 from [Patala 2013] 
        % Inputs:
        %   quat        Disorientation quaternion
        %   region      Region of the FZ (from disorientationregion())
        %   forceaxes   true -> assign the axes even for low symmetry configurations
        % Outputs:
        %   sym    Symmetry group
        %   thrange     Angle of the fundamental zone
        %   x,y,z       Miller indices of the cartesian symmetry axes
        %   x2          Miller indices of the second symmetry axis at angle thrange from x
        function [sym,thrange,x,y,z,x2] = boundaryplanesymmetry(quat,region,forceaxes)
            q0 = quat(1);
            q1 = quat(2);
            q2 = quat(3);
            q3 = quat(4);            

            % Symmetry axes
            if strcmp(region,'OAB')         % Surfaces
                sym = 'C2h';
                thrange = pi;
                x = nan(1,3);
                y = nan(1,3);
                z = [q2,-q1,q0];
            elseif strcmp(region,'OBCE')
                sym = 'C2h';
                thrange = pi;
                x = nan(1,3);
                y = nan(1,3);
                z = [(q0+q3)/sqrt(2),(q3-q0)/sqrt(2),-sqrt(2)*q1];
            elseif strcmp(region,'OADE')
                sym = 'C2h';
                thrange = pi;
                x = nan(1,3);
                y = nan(1,3);
                z = [-sqrt(2)*q2,(q0+q1)/sqrt(2),(q1-q0)/sqrt(2)];
            elseif strcmp(region,'CDE')
                sym = 'C2h';
                thrange = pi;
                x = nan(1,3);
                y = nan(1,3);
                z = [q0-q3,q0-q1,q0-q2];

            elseif strcmp(region,'OB')      % Lines
                sym = 'D2h';
                thrange = pi/2;
                x = [q1,-q1,q0];
                y = [q0/sqrt(2),-q0/sqrt(2),-sqrt(2)*q1];        % ERROR in ay from [Patala 2013]
                z = [1/sqrt(2),1/sqrt(2),0];
            elseif strcmp(region,'CE')
                sym = 'D2h';
                thrange = pi/2;
                x = [2*q1,q1+q3,q1+q3];
                y = sqrt(2)*[-(q1+q3),q1,q1];
                z = [0,-1/sqrt(2),1/sqrt(2)];
            elseif strcmp(region,'DE')
                sym = 'D2h';
                thrange = pi/2;
                x = [q1+q2,2*q2,q1+q2];
                y = sqrt(2)*[q2,-(q1+q2),q2];
                z = [1/sqrt(2),0,-1/sqrt(2)];
            elseif strcmp(region,'OE')
                sym = 'D3d';
                thrange = pi/3;
%                 x = 1/sqrt(2)*[q0+q1,q1-q0,-2*q1];
%                 y = 1/sqrt(6)*[-3*q1+q0,3*q1+q0,-2*q0];
                x = 1/sqrt(6)*[-3*q1+q0,3*q1+q0,-2*q0];     % Swapping x=ay, y=-ax from [Patala 2013]
                y = -1/sqrt(2)*[q0+q1,q1-q0,-2*q1];
                z = 1/sqrt(3)*[1,1,1];
            elseif strcmp(region,'OA')
                sym = 'D4h';
                thrange = pi/4;
                x = [0,q0,q1];
                y = [0,-q1,q0];
                z = [1,0,0];
            elseif strcmp(region,'*AC')
                sym = 'C2h';
                thrange = pi;
                x = nan(1,3);
                y = nan(1,3);
                z = [0,1/sqrt(2),1/sqrt(2)];
            elseif strcmp(region,'AB')  % = OAB
                sym = 'C2h';
                thrange = pi;
                x = nan(1,3);
                y = nan(1,3);
                z = [q2,-q1,q0];
            elseif strcmp(region,'AD')  % = OADE
                sym = 'C2h';
                thrange = pi;
                x = nan(1,3);
                y = nan(1,3);
                z = [-sqrt(2)*q2,(q0+q1)/sqrt(2),(q1-q0)/sqrt(2)];
            elseif strcmp(region,'BC')  % = OBCE
                sym = 'C2h';
                thrange = pi;
                x = nan(1,3);
                y = nan(1,3);
                z = [(q0+q3)/sqrt(2),(q3-q0)/sqrt(2),-sqrt(2)*q1];
            elseif strcmp(region,'CD')  % = CDE
                sym = 'C2h';
                thrange = pi;
                x = nan(1,3);
                y = nan(1,3);
                z = [q0-q3,q0-q1,q0-q2];

            elseif strcmp(region,'B')       % Points
                sym = 'D2h';
                thrange = pi/2;
                x = [q1,-q1,q0];
                y = [q0/sqrt(2),-q0/sqrt(2),-sqrt(2)*q1];        % ERROR in ay from [Patala, 2013]
                z = [1/sqrt(2),1/sqrt(2),0];
            elseif strcmp(region,'E')
                sym = 'D6h';
                thrange = pi/6;
                x = 1/sqrt(6)*[2,-1,-1];
                y = [0,1/sqrt(2),-1/sqrt(2)];
                z = 1/sqrt(3)*[1,1,1];
            elseif strcmp(region,'A')
                sym = 'D8h';
                thrange = pi/8;
                x = [0,cos(pi/8),sin(pi/8)];
                y = [0,-sin(pi/8),cos(pi/8)];
                z = [1,0,0];
            elseif strcmp(region,'C')
                sym = 'D4h';
                thrange = pi/4;
                x = [1,0,0];
                y = [0,1/sqrt(2),1/sqrt(2)];
                z = [0,-1/sqrt(2),1/sqrt(2)];
            elseif strcmp(region,'D')           % Equivalent to point B
                sym = 'D2h';
                thrange = pi/2;
                x = [-q1,q0,-q1];
                y = [q0/sqrt(2),sqrt(2)*q1,q0/sqrt(2)];
                z = [1/sqrt(2),0,-1/sqrt(2)];
            elseif strcmp(region,'O')
                sym = 'Oh';
                thrange = 0;        % Arbitrarily set to zero
                x = [1,0,0];
                y = [0,1,0];
                z = [0,0,1];
            
            else        % Interior or ABCD(~AC)
                sym = 'Ci';
                thrange = 2*pi;
                [x,y,z] = deal(nan(1,3));
            end

            % Setting [x,y,z] even if they are not symmetry axes
            if isnan(z(1))
                z = [q1,q2,q3]/norm([q1,q2,q3]);
            end

            % Calculating the second axis of symmetry, x2
            if ~(thrange==pi || thrange==2*pi)   % Defining the second symmetry axis
                x2 = [-((-y(3)*z(2)*cos(thrange)+y(2)*z(3)*cos(thrange)+x(3)*z(2)*sin(thrange)-x(2)*z(3)*sin(thrange))/...
                         (x(3)*y(2)*z(1)-x(2)*y(3)*z(1)-x(3)*y(1)*z(2)+x(1)*y(3)*z(2)+x(2)*y(1)*z(3)-x(1)*y(2)*z(3))),...
                        -((y(3)*z(1)*cos(thrange)-y(1)*z(3)*cos(thrange)-x(3)*z(1)*sin(thrange)+x(1)*z(3)*sin(thrange))/...
                         (x(3)*y(2)*z(1)-x(2)*y(3)*z(1)-x(3)*y(1)*z(2)+x(1)*y(3)*z(2)+x(2)*y(1)*z(3)-x(1)*y(2)*z(3))),...
                        -((y(2)*z(1)*cos(thrange)-y(1)*z(2)*cos(thrange)-x(2)*z(1)*sin(thrange)+x(1)*z(2)*sin(thrange))/...
                        (-x(3)*y(2)*z(1)+x(2)*y(3)*z(1)+x(3)*y(1)*z(2)-x(1)*y(3)*z(2)-x(2)*y(1)*z(3)+x(1)*y(2)*z(3)))];
            else
                x2 = nan(1,3);
            end
            
            % Force the assignment of axes even for low symmetry configurations
            if forceaxes    
                if isnan(x(1))                                  
                    if z(3) < 1         % Two vectors perpendicular to nboundary
                        y = [ z(2)/sqrt(z(1)^2+z(2)^2),...
                                -z(1)/sqrt(z(1)^2+z(2)^2),...
                                 0];
                    else
                        y = [1,0,0];
                    end
                    x = cross(y,z);
                    x2 = -x;
                end
            end
            
            % Checking for errors
            tol = 20*eps;
            if abs(dot(x,y))>tol || abs(dot(x,z))>tol || abs(dot(y,z))>tol || abs(dot(z,x2))>tol || ...
               abs(1-sqrt(sum(x.^2)))>tol || abs(1-sqrt(sum(y.^2)))>tol || abs(1-sqrt(sum(z.^2)))>tol || abs(1-sqrt(sum(x2.^2)))>tol
                error('Error with the generated axes');
            end
        end



           
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Symmetry axes - main eigenvector
        % To determine the symmetry axes of a uniaxial load (or stress main axis) 
        % Inputs:
        %   nboundary       Miller indices of the boundary plane
        %   fz              Region of the disorientation FZ (from disorientationregion())
        %   symdichromatic  Point group of the boundary plane
        %   fzboundary      Region of the boundary plane FZ (from boundaryplaneregion())
        %   thrangeboundary theta range of the boundary plane FZ
        %   symx,symy,symz,symx2    Symmetry axes of the boundary plane FZ (from boundaryplanessymmetry())
        %   forceaxes       true -> assign the axes even for low symmetry configurations
        % Outputs:
        %   sym             Point group of the loading axis
        %   thrange         Angle of the fundamental zone
        %   x,y,z           Miller indices of the cartesian symmetry axes
        %   x2              Miller indices of the second symmetry axis at angle thrange from x
        function [sym,thrange,x,y,z,x2] = stressaxissymmetry(nboundary,fz,symdichromatic,fzboundary,thrangeboundary,symx,symy,symz,symx2,forceaxes)
            
            % Symmetry axes
            if strcmp(fz,'O')                               % Single crystal
                sym = 'Oh';
                thrange = 0;
                z = nboundary;
                y = nan(1,3);
                x = nan(1,3);
            elseif strcmp(symdichromatic,'Ci')                 % Boundary plane with Ci symmetry
                sym = 'Ci';
                thrange = 2*pi;
                [x,y,z] = deal(nan(1,3));
            elseif strcmp(symdichromatic,'C2h')                % Boundary plane with C2h symmetry   
                if strcmp(fzboundary,'O')               % Origin
                    sym = 'C2h';
                    thrange = thrangeboundary;
                    z = symz;
                    y = symy;
                    x = symx;
                elseif strcmp(fzboundary,'AB') ||  ...  % Along the mirror plane of C2h
                       strcmp(fzboundary,'A') || strcmp(fzboundary,'B')   
                    sym = 'C2h*';
                    thrange = pi;
                    z = nboundary;
                    y = symz;
                    x = cross(y,z);
                else                            % Interior
                    sym = 'Ci';
                    thrange = 2*pi;
                    [x,y,z] = deal(nan(1,3));
                end


            elseif strcmp(symdichromatic,'D3d')                % Boundary plane with D3d symmetry, i.e. OE   
                if strcmp(fzboundary,'O')         % Vertices
                    sym = symdichromatic;
                    thrange = thrangeboundary;
                    z = symz;
                    y = symy;
                    x = symx;
                elseif strcmp(fzboundary,'A')
                    sym = 'C2h*';
                    thrange = pi;
                    z = symx;
                    y = symy;
                    x = -symz;
                elseif strcmp(fzboundary,'B')
                    sym = 'C2h*';
                    thrange = pi;
                    z = symx2;
                    y = cross(symz,symx2);
                    x = -symz;
                elseif strcmp(fzboundary,'OA')    % Edges
                    sym = 'C2h*';
                    thrange = pi;
                    z = nboundary;
                    y = symy;
                    x = cross(y,z);
                elseif strcmp(fzboundary,'OB')
                    sym = 'C2h*';
                    thrange = pi;
                    z = nboundary;
                    y = cross(symz,symx2);
                    x = cross(y,z);
                elseif strcmp(fzboundary,'AB')
                    sym = 'Ci';
                    thrange = 2*pi;
                    z = nboundary;
                    y = symz;
                    x = cross(y,z);
                else                                % Interior
                    sym = 'Ci';
                    thrange = 2*pi;
                    [x,y,z] = deal(nan(1,3));
                end
            else                                            % Boundary planes with higher symmetry and inversion symmetry
                if strcmp(fzboundary,'O')         % Vertices
                    sym = symdichromatic;
                    thrange = thrangeboundary;
                    z = symz;
                    y = symy;
                    x = symx;
                elseif strcmp(fzboundary,'A')
                    sym = 'D2h';
                    thrange = pi/2;
                    z = symx;
                    y = symy;
                    x = -symz;
                elseif strcmp(fzboundary,'B')
                    sym = 'D2h';
                    thrange = pi/2;
                    z = symx2;
                    y = cross(symz,symx2);
                    x = -symz;
                elseif strcmp(fzboundary,'OA')    % Edges
                    sym = 'C2h*';
                    thrange = pi;
                    z = nboundary;
                    y = symy;
                    x = cross(y,z);
                elseif strcmp(fzboundary,'OB')
                    sym = 'C2h*';
                    thrange = pi;
                    z = nboundary;
                    y = cross(symz,symx2);
                    x = cross(y,z);
                elseif strcmp(fzboundary,'AB')
                    sym = 'C2h*';
                    thrange = pi;
                    z = nboundary;
                    y = symz;
                    x = cross(y,z);
                else                                % Interior
                    sym = 'Ci';
                    thrange = 2*pi;
                    [x,y,z] = deal(nan(1,3));
                end
            end

            % Force the assignment of axes even for low symmetry configurations
            if forceaxes
                if isnan(z(1))
                    z = nboundary;
                end
                if isnan(x(1))
                    if all(abs(symz-nboundary)<=10*eps)  % symz == nboundary
                        y = symy;
                        x = symx;
                    else
                        y = cross(symz,nboundary);
                        y = y/norm(y);
                        x = cross(y,z);
                    end
                end
            end

            % Defining the second symmetry axis, x2
            if ~(thrange==0 || thrange==2*pi)   
                x2 = [-((-y(3)*z(2)*cos(thrange)+y(2)*z(3)*cos(thrange)+x(3)*z(2)*sin(thrange)-x(2)*z(3)*sin(thrange))/...
                         (x(3)*y(2)*z(1)-x(2)*y(3)*z(1)-x(3)*y(1)*z(2)+x(1)*y(3)*z(2)+x(2)*y(1)*z(3)-x(1)*y(2)*z(3))),...
                        -((y(3)*z(1)*cos(thrange)-y(1)*z(3)*cos(thrange)-x(3)*z(1)*sin(thrange)+x(1)*z(3)*sin(thrange))/...
                         (x(3)*y(2)*z(1)-x(2)*y(3)*z(1)-x(3)*y(1)*z(2)+x(1)*y(3)*z(2)+x(2)*y(1)*z(3)-x(1)*y(2)*z(3))),...
                        -((y(2)*z(1)*cos(thrange)-y(1)*z(2)*cos(thrange)-x(2)*z(1)*sin(thrange)+x(1)*z(2)*sin(thrange))/...
                        (-x(3)*y(2)*z(1)+x(2)*y(3)*z(1)+x(3)*y(1)*z(2)-x(1)*y(3)*z(2)-x(2)*y(1)*z(3)+x(1)*y(2)*z(3)))];
            else
                x2 = nan(1,3);
            end

            % Checking for errors
            tol = 20*eps;
            if abs(dot(x,y))>tol || abs(dot(x,z))>tol || abs(dot(y,z))>tol || abs(dot(z,x2))>tol || ...
               abs(1-sqrt(sum(x.^2)))>tol || abs(1-sqrt(sum(y.^2)))>tol || abs(1-sqrt(sum(z.^2)))>tol || abs(1-sqrt(sum(x2.^2)))>tol || ...
               abs(dot(z,nboundary)-1)>tol
                error('Error with the generated axes');
            end
        end
    end
end
