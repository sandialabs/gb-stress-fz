%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to map a grain boundary misorientation into the equivalent
% disorientation in the fundamental zone.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Inputs
% euleranglesI      List of euler angles (ZXZ) for each grain
% euleranglesII
%
%%% Outputs:
% dis               Disorientation structure:
%   refgrain        Grain taken as the reference I (1 or 2)
%   fzdichromatic   Location in the disorientation fundamental zone
%   symdichromatic  Point group of the bicrystal
%   w               Weights (non applicable)
%   euleranglesI    Euler angles of grain I [0,0,0]
%   euleranglesII   Euler angles of grain II (ZXZ)
%   l               Disorientation axis
%   th              Disorientation angle [rad]
%   rv              Rodrigues vector
%   symx,symy,symz,symx2    Axes in the symmetry coordinate system
%   thrangeboundary Angle range of the boundary plane FZ
% r02disfz          Rotation matrix from the original to the disorientation coordinate system in the FZ
%

function [dis,r02disfz] = toFZ_disorientation(euleranglesI,euleranglesII)
    assert(size(euleranglesI,1) == size(euleranglesII,1) && ...
           size(euleranglesI ,2) == 3 && ...
           size(euleranglesII,2) == 3,'euleranglesI and euleranglesII must be nx3 matrices');

    % Preallocating
    n_in = size(euleranglesI,1);                                    % Number of points

    [euleranglesI_FZ,euleranglesII_FZ,rv_FZ,l_FZ] = deal(nan(n_in,3));
    [th_FZ,refgrain] = deal(nan(n_in,1));
    [r02disfz,fzdichromatic,rdis] = deal(cell(n_in,1));
    
    %%%%% Finding the FZ disorientation
    symrot = symmetryrotations('Oh');                               % Proper symmetry rotations
    nsym = length(symrot);                                          % Number of proper symmetry rotations
    for n = 1:n_in

        % 0 coordinate system - initial rotation matrices
        rI0  = trans.rEuler( euleranglesI(n,1), euleranglesI(n,2), euleranglesI(n,3))';
        rII0 = trans.rEuler(euleranglesII(n,1),euleranglesII(n,2),euleranglesII(n,3))';

        % Misorientations referenced to each grain
        MI  = rI0'*rII0;                                               
        MII = rII0'*rI0;
        % Note: transformations are for misorientations as active rotations;
        % alternatives if these were passive rotations instead are:
        % MI  = rII0*rI0';            % rI0*(rI0'*rII0)*rI0';                                           
        % MII = rI0*rII0';            % rII0*(rII0'*rI0)*rII0';


        % Testing symmetry rotations
        l = nan(2*nsym^2,3);                                        % Initializing
        th = nan(2*nsym^2,1);
        c = 1;                                                      % Counter
        
        % ---> Taking grain I as the reference
        for i = 1:nsym
            for j = 1:nsym
                % Symmetric misorientations
                Ms = symrot{i}*MI*symrot{j}';
    
                % Axis and angle representation
                [l(c,:),th(c)] = trans.r2axisangle(Ms);
                c = c+1;
            end
        end
        % ---> Taking grain II as the reference
        for i = 1:nsym
            for j = 1:nsym
                % Symmetric misorientations
                Ms = symrot{i}*MII*symrot{j}';
    
                % Axis and angle representation
                [l(c,:),th(c)] = trans.r2axisangle(Ms);
                c = c+1;
            end
        end
        % Conditions with the lowest th (i.e. disorientations)
        [th_sorted,id_sorted] = sort(th,'ascend');
        dis_ids = id_sorted(th_sorted-th_sorted(1) <= 40*eps);      % Large tolerance was found to be required    

        % Finding the disorientation within the FZ
        for i = 1:length(dis_ids)
            idi = dis_ids(i);
            rv = trans.axisangle2rodrigues(l(idi,:),th(idi));       % Rodrigues vector

            tol = 2*eps;                                            % Tolerance
            if all(rv >= -tol) && ...                               % FZ conditions
               rv(2) >= rv(3)-tol && ...
               rv(1) >= rv(2)-tol             

                % Region of the FZ
                fzn = regionSO3('Oh','Oh',rv,1);
                if ~(strcmp(fzn,'D') || ...                         % Ignore these, as point D = B
                     strcmp(fzn,'AD') || ...
                     strcmp(fzn,'CD'))
                    
                    fzdichromatic{n} = fzn{1};                      % Assign FZ
    
                    % Fundamental zone parameters
                    rv_FZ(n,:) = rv;
                    l_FZ(n,:) = l(idi,:);
                    th_FZ(n) = th(idi);
        
                    % Recovering the rotations angles
                    rdis{n} = trans.raxisangle(l(idi,:),th(idi));   % Disorientation
                    euleranglesI_FZ(n,:) = [0,0,0];                 % Euler angles
                    [euleranglesII_FZ(n,1),euleranglesII_FZ(n,2),euleranglesII_FZ(n,3)] = trans.r2euler(rdis{n});
                    
                    % Rotation from 0 to dis coordinate systems
                    if idi <= nsym^2                                % Grain I as reference
                        symi = symrot{floor((idi-1)/nsym)+1};       % Symmetry operation i
                        symj = symrot{mod(idi-1,nsym)+1};           % Symmetry operation j
                        r02disfz{n} = symj * rII0';                 % Rotation matrix for GBplane and stress
                        refgrain(n) = 1;                            % Reference grain
                    else                                            % Grain II as reference
                        symi = symrot{floor((idi-nsym^2-1)/nsym)+1};
                        symj = symrot{mod(idi-nsym^2-1,nsym)+1};
                        r02disfz{n} = symj * rI0';
                        refgrain(n) = 2;
                    end

                    % Exit loop
                    break
                end
            end
        end 
    end
    
    % Symmetry and symmetry axes for each disorientation - Table 21 of [Patala 2013]
    symdichromatic = cell(size(rv_FZ,1),1);                 % Symmetry
    thrangeboundary = nan(size(rv_FZ,1),1);                 % th range for the boundary FZ
    [symx,symy,symz,symx2] = deal(nan(size(rv_FZ,1),3));    % Symmetry axes
    for i = 1:size(rv_FZ,1)
        quat = trans.rodrigues2quaternion(rv_FZ(i,:));      % Quaternion
        [symdichromatic{i},thrangeboundary(i),symx(i,:),symy(i,:),symz(i,:),symx2(i,:)] = fccFZ.boundaryplanesymmetry(quat,fzdichromatic{i},true);
    end


    %%%%% Building dis structure
    dis.refgrain = refgrain;
    dis.fzdichromatic = fzdichromatic;
    dis.symdichromatic = symdichromatic;
    dis.w = ones(size(fzdichromatic));
    dis.euleranglesI  = euleranglesI_FZ ;
    dis.euleranglesII = euleranglesII_FZ;
    dis.l = l_FZ;
    dis.th = th_FZ;
    dis.rv = rv_FZ;
    dis.symx = symx;
    dis.symy = symy;
    dis.symz = symz;
    dis.symx2 = symx2;
    dis.thrangeboundary = thrangeboundary;

end