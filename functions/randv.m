%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions to generate random vectors and stress tensors.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Functions:
% randEuler
% randVector
% randStress
%

classdef randv      
    methods(Static)
        
        %%%%% Random outputs
        % Random Euler angles (ZXZ) for a homogeneous orientation distribution
        function eulerangles = randEuler(n)
            assert(all(size(n)==[1,1]) && mod(n,1)==0, 'n must be a single integer');
            eulerangles = nan(n,3);
            for i = 1:n
                quat = compact(randrot(1,1));               % Random quaternion
                r = trans.quat2r(quat);                     % Rotation matrix
                [phi1,Phi,phi2] = trans.r2euler(r);         % Euler angles
                eulerangles(i,:) = [phi1,Phi,phi2];
            end
        end


        % Random unit vectors for a homogeneous orientation distribution
        % Output modes:
        %   vector          Uniaxial vector
        %   vectorangles    Angles [phi,th] of the uniaxial vector
        function result = randVector(n,mode)
            assert(all(size(n)==[1,1]) && mod(n,1)==0, 'n must be a single integer');
            assert(strcmp(mode,'vector') || ...
                   strcmp(mode,'vectorangles'), "mode must be 'vector' or 'vectorangles'")
            
            vector = normalize(randn(3,n),'norm')';         % Unit vectors
            if strcmp(mode,'vector')                        % Output vectors
                result = vector;
            elseif strcmp(mode,'vectorangles')              % Output angles
                result = nan(n,2);
                for i = 1:n
                    [~,result(i,1),result(i,2)] = trans.cart2sph(vector(i,:));  % Converting to [phi,theta]
                end
            end
        end

        % Random stress tensor in Voigt notation
        % Output modes:
        %   uniaxial        Uniaxial stress (e3 = 1, e1 = e2 = 0)
        %   biaxial         Biaxial stresses (e3 = e2, e1 = 0)
        %   axisymetric     Biaxial stresses (e3 = e2 ~= e1 )
        %   triaxial        Triaxial stress (e3 - e1 = 1) 
        % * [e1,e2,e3] are the eigenvalues such that e1 >= e2 >= e3
        function result = randStress(n,mode)
            assert(strcmp(mode,'uniaxial') || ...
                   strcmp(mode,'biaxial') || ...
                   strcmp(mode,'axisymmetric') || ...
                   strcmp(mode,'triaxial'), "mode must be 'uniaxial', 'biaxial', 'axisymmetric' or 'triaxial'")
            if strcmp(mode,'uniaxial') 
                vector = normalize(randn(3,n),'norm')';         % Unit vectors
                result = nan(n,6);
                for i = 1:n
                    stress = vector(i,:)' * vector(i,:);        % Stress tensor
                    result(i,:) = trans.N2Vstress(stress);
                end
            elseif strcmp(mode,'biaxial') 
                vector = normalize(randn(3,n),'norm')';         % Unit vectors
                stress0 = [0,0,0;0,1,0;0,0,1];
                result = nan(n,6);
                for i = 1:n
                    v1 = vector(i,:)';
                    if abs(v1(3)) < 0.9 % check if v1 is not too close to the z-axis
                        temp_v = [0; 0; 1];
                    else
                        temp_v = [1; 0; 0];
                    end       
                    v2 = normalize(cross(temp_v, v1),'norm');
                    v3 = cross(v1, v2);
                    Q = [v1, v2, v3];                           % Rotation matrix
                    stressR = Q' * stress0 * Q;                 % Rotated stress tensor
                    result(i,:) = trans.N2Vstress(stressR);
                end
            elseif strcmp(mode,'axisymmetric') 
                vector = normalize(randn(3,n),'norm')';         % Unit vectors
                result = nan(n,6);
                for i = 1:n
                    estress = randn(1,2);                       % Eigenstresses
                    estressnorm = estress/abs(diff(estress));   % Normalized
                    if max(abs(estressnorm)) > 1                % Bring back to mean(e)~0
                        estressnorm = estressnorm - max(estressnorm) + 1;
                    end
                    stress0 = [estressnorm(1),0,0;...           % Axisymmetric stress
                               0,estressnorm(2),0;...
                               0,0,estressnorm(2)];
                    v1 = vector(i,:)';
                    if abs(v1(3)) < 0.9 % check if v1 is not too close to the z-axis
                        temp_v = [0; 0; 1];
                    else
                        temp_v = [1; 0; 0];
                    end       
                    v2 = normalize(cross(temp_v, v1),'norm');
                    v3 = cross(v1, v2);
                    Q = [v1, v2, v3];                           % Rotation matrix
                    stressR = Q' * stress0 * Q;                 % Rotated stress tensor
                    result(i,:) = trans.N2Vstress(stressR);
                end 
            elseif strcmp(mode,'triaxial') 
                result = nan(n,6);
                for i = 1:n
                    stress = rand(1,6);                         % Random stress
                    ev = eig(trans.V2Nstress(stress));          % Eigenvalues
                    result(i,:) = stress/(ev(3)-ev(1));         % Normalized stress
                end
            end
        end

    end
end