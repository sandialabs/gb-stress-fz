%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiple transformations for vectors and tensors.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% General notation:
% r     rotation matrix
% R     rotate an object (vector, tensor)
% _V    Voigt notation
% _N    Normal (full matrix) notation
%
%%% Functions:
% normalize,            
% N2Vstrain,            V2Nstrain,              N2Vstress,          V2Nstress
% sph2cart,             cart2sph,               cart2sphM,          pol2sphM
% rx,                   ry,                     rz,                 
% rVN,                  rVM,
% raxisangle,           rEuler,                 r2euler,            r2axisangle,
% r2rotationvector,     r2rodrigues,            r2quat,             axisangle2rodrigues,
% rodrigues2axisangle,  rodrigues2quaternion,   rodrigues2r,        quat2r,
% vector2skewsymmatrix, skewsymmatrix2vector,
% R,
% RVcompliance,         RVstiffness,            RVstrain,           RVstress,
% direction2Miller,     latexNegativeMiller 
%

classdef trans      
    methods(Static)
        
        %%%%% Basic functions
        function v = normalize(v)
            v = v/norm(v);
        end

        %%%%% Normal vs Voigt notations
        function eV = N2Vstrain(eN)
            eV = [eN(1,1),eN(2,2),eN(3,3),2*eN(2,3),2*eN(1,3),2*eN(1,2)]';
        end

        function eN = V2Nstrain(eV)
            eN = [    eV(1),1/2*eV(6),1/2*eV(5);...
                  1/2*eV(6),    eV(2),1/2*eV(4);...
                  1/2*eV(5),1/2*eV(4),    eV(3)];
        end
        
        function sV = N2Vstress(sN)
            sV = [sN(1,1),sN(2,2),sN(3,3),sN(2,3),sN(1,3),sN(1,2)]';
        end

        function sN = V2Nstress(sV)
            sN = [sV(1),sV(6),sV(5);...
                  sV(6),sV(2),sV(4);...
                  sV(5),sV(4),sV(3)];
        end
       
        %%%%% Vector transformations
        % Spherical to cartesian
        function result = sph2cart(rho,theta,phi)
            result = [rho*sin(phi)*cos(theta);...
                      rho*sin(phi)*sin(theta);...
                      rho*cos(phi)];
        end

        % Cartesian to spherical
        function [rho,theta,phi] = cart2sph(cartvector)
            assert(all(size(cartvector) == [1,3]) || ...
                   all(size(cartvector) == [3,1]),'cartvector must be a 3x1 or 1x3 matrix')
            x = cartvector(1);
            y = cartvector(2);
            z = cartvector(3);
            if x == 0 && y == 0 && z >= 0    % Cases
                rho = z;
                theta = 0;
                phi = 0;
            elseif x == 0 && y == 0 && z < 0 
                rho = abs(z);
                theta = 0;
                phi = pi;
            else                            % General case
                rho = sqrt(x^2+y^2+z^2);
                theta = atan2(y,x);
                phi = acos(z/rho);
            end
        end

        % Cartesian to spherical (series of vectors)
        % cartmat       [x,y,z]
        function [rho,theta,phi] = cart2sphM(cartmat)
            assert(size(cartmat,2) == 3,'cartmat must have 3 columns')
            x = cartmat(:,1);
            y = cartmat(:,2);
            z = cartmat(:,3);

            rho = sqrt(x.^2+y.^2+z.^2);             % General case
            theta = atan2(y,x);
            phi = acos(z./rho);

            icase1 = (x == 0 & y == 0 & z >= 0);  % Case 1
            rho(icase1) = z(icase1);
            theta(icase1) = 0;
            phi(icase1) = 0;
            
            icase2 = (x == 0 & y == 0 & z < 0);   % Case 2 
            rho(icase2) = abs(z(icase2));
            theta(icase2) = 0;
            phi(icase2) = pi;
        end

        % Polar (3D) to spherical
        % polarmat      [r,theta,z]
        function [rho,theta,phi] = pol2sphM(polarmat)
            assert(size(polarmat,2) == 3,'polarmat must have 3 columns')
            r = polarmat(:,1);
            th = polarmat(:,2);
            z = polarmat(:,3);

            x = r.*cos(th);
            y = r.*sin(th);

            [rho,theta,phi] = trans.cart2sphM([x,y,z]);
        end

        %%%%% Rotation matrices (r_)
        % Rotation matrices along each axis
        function result = rx(th)
            result = [       1,        0,        0;... 
                             0,  cos(th),  sin(th);...
                             0, -sin(th),  cos(th)];
        end

        function result = ry(th)
            result = [ cos(th),        0, -sin(th);...
                             0,        1,        0;...
                       sin(th),        0,  cos(th)];
        end

        function result = rz(th)
            result = [ cos(th),  sin(th),        0;...
                      -sin(th),  cos(th),        0;...
                             0,        0,        1];
        end
        
        % Rotations in Voigt notation
        function result = rVN(r)
            assert(all(size(r) == [3,3]), 'r must be a 3x3 matrix');
            result = [r(1,1)^2, r(1,2)^2, r(1,3)^2, r(1,2)*r(1,3), r(1,3)*r(1,1), r(1,1)*r(1,2);...            % Voigt compliance rotation matrix for a rotation matrix r [Auld -  Acoustic fields and waves in solids Volume I]
                      r(2,1)^2, r(2,2)^2, r(2,3)^2, r(2,2)*r(2,3), r(2,3)*r(2,1), r(2,1)*r(2,2);...
                      r(3,1)^2, r(3,2)^2, r(3,3)^2, r(3,2)*r(3,3), r(3,3)*r(3,1), r(3,1)*r(3,2);...
                      2*r(2,1)*r(3,1), 2*r(2,2)*r(3,2), 2*r(2,3)*r(3,3), r(2,2)*r(3,3)+r(2,3)*r(3,2), r(2,1)*r(3,3)+r(2,3)*r(3,1), r(2,2)*r(3,1)+r(2,1)*r(3,2);...
                      2*r(3,1)*r(1,1), 2*r(3,2)*r(1,2), 2*r(3,3)*r(1,3), r(1,2)*r(3,3)+r(1,3)*r(3,2), r(1,3)*r(3,1)+r(1,1)*r(3,3), r(1,1)*r(3,2)+r(1,2)*r(3,1);...
                      2*r(1,1)*r(2,1), 2*r(1,2)*r(2,2), 2*r(1,3)*r(2,3), r(1,2)*r(2,3)+r(1,3)*r(2,2), r(1,3)*r(2,1)+r(1,1)*r(2,3), r(1,1)*r(2,2)+r(1,2)*r(2,1)];
        end
        
        function result = rVM(r)
            assert(all(size(r) == [3,3]), 'r must be a 3x3 matrix');
            result = [r(1,1)^2, r(1,2)^2, r(1,3)^2, 2*r(1,2)*r(1,3), 2*r(1,3)*r(1,1), 2*r(1,1)*r(1,2);...      % Voigt stiffness rotation matrix for a rotation matrix r [Auld -  Acoustic fields and waves in solids Volume I]  
                      r(2,1)^2, r(2,2)^2, r(2,3)^2, 2*r(2,2)*r(2,3), 2*r(2,3)*r(2,1), 2*r(2,1)*r(2,2);...
                      r(3,1)^2, r(3,2)^2, r(3,3)^2, 2*r(3,2)*r(3,3), 2*r(3,3)*r(3,1), 2*r(3,1)*r(3,2);...
                      r(2,1)*r(3,1), r(2,2)*r(3,2), r(2,3)*r(3,3), r(2,2)*r(3,3)+r(2,3)*r(3,2), r(2,1)*r(3,3)+r(2,3)*r(3,1), r(2,2)*r(3,1)+r(2,1)*r(3,2);...
                      r(3,1)*r(1,1), r(3,2)*r(1,2), r(3,3)*r(1,3), r(1,2)*r(3,3)+r(1,3)*r(3,2), r(1,3)*r(3,1)+r(1,1)*r(3,3), r(1,1)*r(3,2)+r(1,2)*r(3,1);...
                      r(1,1)*r(2,1), r(1,2)*r(2,2), r(1,3)*r(2,3), r(1,2)*r(2,3)+r(1,3)*r(2,2), r(1,3)*r(2,1)+r(1,1)*r(2,3), r(1,1)*r(2,2)+r(1,2)*r(2,1)];
        end

        %%%%% Rotations - converting between notations

        % Rotation matrix from angle-axis notation
        function result = raxisangle(l,th)
            if th ~= 0
                % l = axis [lx, ly, lz], th = rotation angle
                lx = l(1)/norm(l); ly = l(2)/norm(l); lz = l(3)/norm(l);
                c = cos(th);
                s = -sin(th);           % + or - changes between active/passive rotation
                result = [lx^2*(1-c)+c,     lx*ly*(1-c)-lz*s,   lx*lz*(1-c)+ly*s;...
                          lx*ly*(1-c)+lz*s, ly^2*(1-c)+c,       ly*lz*(1-c)-lx*s;...
                          lx*lz*(1-c)-ly*s, ly*lz*(1-c)+lx*s,   lz^2*(1-c)+c];
            else
                result = [1,0,0;0,1,0;0,0,1];
            end
        end
        
        % Rotation by Euler angles (ZXZ)
        function result = rEuler(phi1,Phi,phi2)
            result = trans.rz(phi2) * trans.rx(Phi) * trans.rz(phi1);
        end

        % Rotation matrix -> Euler angles
        function [phi1, Phi, phi2] = r2euler(R)
            assert(all(size(R) == [3,3]) && abs(1-det(R))<=10*eps, 'Input must be a 3x3 rotation matrix');
            if abs(R(3,3)) < 1 - eps                        % Extract angles
                Phi = acos(R(3,3));
                phi1 = atan2(R(3,1), -R(3,2));
                phi2 = atan2(R(1,3), R(2,3));
            else                                            % Special case: Gimbal lock when R(3,3) is ±1
                Phi = acos(R(3,3));                         % 0 or pi
                phi1 = atan2(R(1,2), R(1,1));
                phi2 = 0;                                   % Can be arbitrarily set to 0
            end
        
            phi1 = mod(phi1, 2*pi);                         % Ensure angles are in the range [0, 2π)
            Phi = mod(Phi, 2*pi);
            phi2 = mod(phi2, 2*pi);
        end

        % Rotation matrix -> Axis-angle
        function [l,th] = r2axisangle(R)
            assert(all(size(R) == [3,3]) && abs(1-det(R))<=10*eps, 'Input must be a 3x3 rotation matrix');
            R = R';                                                 % To adapt it for passive rotations
            th = acos((trace(R) - 1) / 2);
            if abs(th) < eps                                        % Special case: 0 deg
                l = [0,0,0];
            elseif abs(th-pi) < eps                                 % Special case: 180 deg
                [V, D] = eig(R);
                [~,idx] = min(abs(diag(D)-1));
                l = normalize(V(:,idx),'norm')';
            else                                                    % General case
                S = (R - R')/(2 * sin(th));
                l = [S(3,2), S(1,3), S(2,1)];
                l = l/norm(l);
            end
        end

        % Rotation matrix -> Rotation vector ([axis]*angle)
        function rv = r2rotationvector(R)
            assert(all(size(R) == [3,3]) && abs(1-det(R))<=10*eps, 'Input must be a 3x3 rotation matrix');
            R = R';                                                 % To adapt it for passive rotations
            if trace(R)+1<0 && trace(R)+1>-4*eps                    % To fix numerical errorwhen trace(R)~-1
                theta = acos(-1);
            else
                theta = acos((trace(R) - 1) / 2);
            end
            if abs(theta) < eps || all(abs(R-R')<4*eps,'all')       % Special cases: 0 and 360 deg
                rv = [0,0,0];
            else                                                    % General case
                S = (R - R')/(2 * sin(theta));
                axis = [S(3,2), S(1,3), S(2,1)];
                rv = axis*theta;
            end
        end

        % Rotation matrix -> Rodrigues vector
        function rv = r2rodrigues(R)
            assert(all(size(R) == [3,3]) && abs(1-det(R))<=10*eps, 'Input must be a 3x3 rotation matrix');
            R = R';                                                 % To adapt it for passive rotations
            if trace(R)+1<0 && trace(R)+1>-4*eps                    % To fix numerical errorwhen trace(R)~-1
                theta = acos(-1);
            else
                theta = acos((trace(R) - 1) / 2);
            end
            if abs(theta) < eps || all(abs(R-R')<4*eps,'all')       % Special cases: 0 and 360 deg
                rv = [0,0,0];
            else                                                    % General case
                S = (R - R')/(2 * sin(theta));
                axis = [S(3,2), S(1,3), S(2,1)];
                rv = axis*tan(theta/2);
            end
        end

        % Rotation matrix -> quaternion
        function quat = r2quat(R)
            assert(all(size(R) == [3,3]) && abs(1-det(R))<=10*eps, 'Input must be a 3x3 rotation matrix');
            [l,th] = trans.r2axisangle(R');
            quat = [cos(th/2),sin(th/2)*l];
        end

        % Axis-angle -> Rodrigues vector
        function rodriguesvector = axisangle2rodrigues(l,th)
            assert(all(size(l) == [1,3]) && (abs(1-norm(l))<=2*eps || all(l == [0,0,0])), 'l must be a unit 1x3 vector');
            assert(th >= 0 && th < 2*pi, 'th range: 0 <= th < 2*pi');
            if th ~= pi 
                rodriguesvector = l*tan(th/2);
            else
                rodriguesvector = l*inf;
            end
        end

        % Rodrigues vector -> Axis-angle
        function [l,th] = rodrigues2axisangle(rodriguesvector)
            assert(all(size(rodriguesvector) == [1,3]), 'rodriguesvector must be a 1x3 vector');
            mag = norm(rodriguesvector);
            if mag == 0
                l = [0,0,1];
                th = 0;
            elseif mag == inf
                l = [0,0,1];
                th = pi;
            else
                l = rodriguesvector/mag;
                th = 2*atan(mag);
            end
        end

        % Rodrigues vector -> quaternion
        function quat = rodrigues2quaternion(rodriguesvector)
            assert(all(size(rodriguesvector) == [1,3]), 'rodriguesvector must be a 1x3 vector');
            theta = norm(rodriguesvector);
            if theta == 0                                   % Identity rotation
                quat = [1, 0, 0, 0];  
            elseif isinf(theta)                             % 180-degree rotation
                quat = [0, rodriguesvector/norm(rodriguesvector)];  
            else                                            % General case
                w = 1 / sqrt(1 + theta^2);
                v = w * rodriguesvector;
                quat = [w,v];
            end
        end

        % Rodrigues vector -> rotation matrix
        function r = rodrigues2r(rodriguesvector)
            assert(all(size(rodriguesvector) == [1,3]), 'rodriguesvector must be a 1x3 vector');
            [l,th] = trans.rodrigues2axisangle(rodriguesvector);
            r = trans.raxisangle(l,th);
        end

        % Quaternion -> rotation matrix
        function r = quat2r(quat)
            assert(all(size(quat) == [1,4]), 'quat must be a 1x4 vector [w,x,y,z]');
            quat = quat / norm(quat);       
            w = quat(1); 
            x = quat(2); 
            y = quat(3); 
            z = quat(4);
        
            ww = w*w; xx = x*x; yy = y*y; zz = z*z;         % Precompute for efficiency
            wx = w*x; wy = w*y; wz = w*z;
            xy = x*y; xz = x*z; yz = y*z;
        
            r = [ww+xx-yy-zz,   2*(xy-wz), 2*(xz + wy);     % Rotation matrix
                 2*(xy + wz), ww-xx+yy-zz, 2*(yz - wx);
                 2*(xz - wy),   2*(yz+wx), ww-xx-yy+zz];
        end

        %%%%% Transforming vector <--> matrix
        % Vector to skew-symetric matrix
        function M = vector2skewsymmatrix(v)
            assert(all(size(v) == [3,1]), 'v must be a 3x1 vector');
            M = [    0 -v(3)  v(2);
                  v(3)     0 -v(1);
                 -v(2)  v(1)     0];
        end

        % Skew-symetric matrix to vector
        function v = skewsymmatrix2vector(M)
            assert(all(size(M) == [1,3]) && ...
                   issymmetric(M,'skew'), 'M must be a 3x3 skew-symmetric matrix');
            v = [M(3,2),M(1,3),M(2,1)]';
        end


        %%%%% Rotating the actual tensors (R_)
        % Normal notation
        function result = R(tensor,r)                       % Normal notation (for stress and strain)
            assert(all(size(tensor) == [3,3]), 'tensor must be a 3x3 vector');
            assert(all(size(r) == [3,3]), 'r must be a 3x3 matrix');
            result = r * tensor * r';                                             
        end

        % Voigt notation
        function result = RVcompliance(S,r)                 % Compliance Voigt matrix rotation
            assert(all(size(S) == [6,6]), 'S must be a 6x6 vector');
            assert(all(size(r) == [3,3]), 'r must be a 3x3 matrix');
            rot = trans.rVN(r);
            result = rot * S * rot';                                             
        end
        
        function result = RVstiffness(C,r)                  % Stiffness Voigt matrix rotation
            assert(all(size(C) == [6,6]), 'C must be a 6x6 vector');
            assert(all(size(r) == [3,3]), 'r must be a 3x3 matrix');
            rot = trans.rVM(r);
            result = rot * C * rot';                                             
        end

        function result = RVstrain(eV,r)                    % Stress Voigt vector rotation
            assert(all(size(eV) == [6,1]), 'eV must be a 6x1 vector');
            assert(all(size(r) == [3,3]), 'r must be a 3x3 matrix');
            rot = trans.rVN(r);
            result = rot * eV;                                             
        end

        function result = RVstress(sV,r)                    % Stress Voigt vector rotation
            assert(all(size(sV) == [6,1]), 'sV must be a 6x1 vector');
            assert(all(size(r) == [3,3]), 'r must be a 3x3 matrix');
            rot = trans.rVM(r);
            result = rot * sV;                                             
        end

        
        %%%%% Miller indices

        % Unit vector (float) to Miller indices (text)
        % Inputs:
        % v         Vector
        % mode      0 -> numeric output
        %           1 -> string (latex interpreter) output
        % *NOTE, format example: xlabel(['$\mathsf{',trans.direction2Miller(v,1),'}$'],'Interpreter','latex')
        function result = direction2Miller(v,mode)
            assert(all(size(v) == [1,3]), 'v must be a 1x3 vector');

            v = v/max(abs(v));

            % Finding the best Miller indices iterating up to 99
            tolerance = 1e-4;           % Tolerance
            max_index_limit = 99;
            bestC = max_index_limit;    % Default to the largest possible scale factor
        
            smallest_error = 1000;
            for C = 1:max_index_limit
                H_C_float = C*v;
                H_C_int = round(H_C_float);
                errorC = abs(H_C_float - H_C_int);  % Difference between float and integer
                if all(errorC < tolerance)
                    bestC = C;
                    break;              % Stop iteration
                elseif sum(abs(errorC)) < smallest_error    % Best so far
                    smallest_error = sum(abs(errorC));
                    bestC = C;
                end
            end
            M = round(bestC*v);         % Miller indices

            if mode == 0
                result = M;
            elseif mode == 1
                % Converting to string
                if all(abs(M)<10)
                    delimiter = '\,';
                else
                    delimiter = '\,\,';
                end
                result = ['[\,',trans.latexNegativeMiller(M(1)),delimiter,...
                    trans.latexNegativeMiller(M(2)),delimiter,...
                    trans.latexNegativeMiller(M(3)),'\,]'];
            end
        end
    
        % Convert a negative character x to \overline{x} for the latex interpreter
        function result = latexNegativeMiller(x)
            if x < 0 
                result = ['\overline{',num2str(abs(x)),'}'];
            else
                result = num2str(x);
            end
        end


    end
end



