%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quaternion and octonion operations.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Functions:
% Qconj
% Qnorm
% Qmult
% Qdot
% Qangle
% Qrotz
% quat2oct
% Oconj
% Onorm
% Omult
% Odot
% Oangle
% Orotz
%

classdef qo      
    methods(Static)

        %%%%%%%%%%%%%%%%
        %% Quaternion %%
        %%%%%%%%%%%%%%%%

        % Conjugate
        function res = Qconj(p) 
            p0 = p(1);
            pv = p(2:end);
            res = [p0,-pv];
        end

        % Norm
        function res = Qnorm(p) 
            res = norm(p);
        end

        % Multiplication
        function res = Qmult(p,q) 
            p0 = p(1);
            pv = p(2:end);
            q0 = q(1);
            qv = q(2:end);

            res = [p0*q0 - dot(pv,qv),...
                   p0*qv + q0*pv + cross(pv,qv)];
        end

        % Dot product
        function res = Qdot(p,q) 
            res = dot(p,q);
        end

        % Angle between two quaternions
        function res = Qangle(p,q) 
            res = 2*acos(abs(qo.Qdot(p,q)));
        end

        % Rotate a quaternion by angle th along the z-axis
        function res = Qrotz(p,th)
            qrotz = [cos(th/2),0,0,sin(th/2)];
            res = qo.Qmult(qrotz,p);
        end

        %%%%%%%%%%%%%%%%
        %%  Octonion  %%
        %%%%%%%%%%%%%%%%

        % Building an octonion from two quaternions
        function res = quat2oct(p,q) 
            res = 1/sqrt(2) * [p,q];
        end

        % Conjugate
        function res = Oconj(p) 
            p0 = p(1);
            pv = p(2:end);
            res = [p0,-pv];
        end

        % Norm
        function res = Onorm(p) 
            res = norm(p);
        end

        % Multiplication - via Cayley-Dickson construction
        function res = Omult(p,q) 
            a = p(1:4);
            b = p(5:8);
            c = q(1:4);
            d = q(5:8);

            res = [qo.Qmult(a,c) - qo.Qmult(qo.Qconj(d),b),...
                   qo.Qmult(d,a) + qo.Qmult(b,qo.Qconj(c))];
        end

        % Dot product
        function res = Odot(p,q) 
            res = dot(p,q);
        end

        % Angle between two quaternions
        function res = Oangle(p,q) 
            res = 2*acos(abs(qo.Odot(p,q)));
        end

        % Rotate an octonion by angle th along the z-axis
        function res = Orotz(p,th)
            res = [qo.Qrotz(p(1:4),th),qo.Qrotz(p(5:8),th)];
        end

    end
end