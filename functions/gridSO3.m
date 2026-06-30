%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to create a grid in SO(3) orientation space.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:
% pg1, pg2      Point group symmetry of each element
% steps         Grid steps in each dimension (every 90deg rotation)
% ge            Grain exchange symmetry? 0 -> no, 1 -> yes
% reg           Determine point within FZ for each point? 0 -> no, 1 -> yes
%               (embedded within gridSO3 so that the FZ does not need to be
%               recalculated for each point)
%
% Outputs:
% rv            Rodrigues vector components [rx,ry,rz]
% region        Region within the FZ for each point
% w             Weights (i.e. number of symmetric points)
% V             Vertices of the Rodrigues-Frank fundamental zone
% F             Faces of the fundamental zone (vertex connectivity)
%

function [rv,region,w,V,F] = gridSO3(pg1,pg2,steps,ge,reg)
    assert(ischar(pg1) && ischar(pg2),'pg1 and pg2 must be point group symmetry strings')
    assert((~ischar(steps) && length(steps)==1 && all(floor(steps)==steps)) || ...
       strcmp(steps,'unique'), ...
       "steps must be an integer or the string 'unique'")
    assert(ge == 0 || ge == 1,'ge must be 0 or 1')
    assert(reg == 0 || reg == 1,'reg must be 0 or 1')

    pgsorted = sort(string({pg1,pg2}));         % Sort point groups
    pg1 = pgsorted(1);
    pg2 = pgsorted(2);

    % Fundamental zone in Rodrigues-Frank space
    [V,F,thmax] = fzSO3(pg1,pg2);             % Vertices and faces      
    
    %%%%%%%%%%%%%%%%% GRID
    if isnumeric(steps)
    
        %%%%% Oh-Oh
        if (ge == 1) && strcmp(pg1,'Oh') && strcmp(pg2,'Oh')
            steps_gIIrotationaxis = steps;
            steps_gIIrotationangle = steps;
               
            % Grain II - rotation through an axis l and disorientation angle th
            lsize = 0;
            for i = 1:steps_gIIrotationangle    % Calculating the size of l
                stepsli = ceil(steps_gIIrotationaxis * i/steps_gIIrotationangle);
                lsize = lsize + 1/2*(stepsli+1)*(stepsli+2);
            end
            
            rv = nan(lsize,3);
            rv001max = tan(pi/8);       % Planes of FZ in Rodrigues vectors
            c = 1;
            for i = 1:steps_gIIrotationangle        
                stepsli = ceil(steps_gIIrotationaxis * i/steps_gIIrotationangle);  % Steps in rotation axis: More steps for larger disorientation angles
                ljv = 0:1/stepsli:1;
                for j = 1:length(ljv)
                    lkv = 0:1/stepsli:ljv(j);
                    for k = 1:length(lkv)
                        % Disorientation axis
                        l = trans.normalize([1,ljv(j),lkv(k)]);
            
                        % Disorientation in Rodrigues space
                        li = l/l(1);                   
                        rv1max = min(li(1)/sum(li),rv001max);                       % Largest possible Rodrigues vector = li*R1max 
                        thmax = 2*atan(norm(li*rv1max));                            % Largest possible misorientation angle for such axis
    
                        th = thmax * i/steps_gIIrotationangle;                      % Disorientation angle
                        rv(c,:) = li/norm(li)*tan(th/2);                            % Rodrigues vector
            
                        % Advancing the counter    
                        c = c+1;
                    end
                end
            end
    
        
        %%%%% D2h-[D4h,D6h,D8h,D3d] - Polygon extruded in z
        elseif (ge == 0) && ...
               ((strcmp(pg1,'D2h') && (strcmp(pg2,'D8h') || strcmp(pg2,'D6h') || strcmp(pg2,'D4h') || strcmp(pg2,'D3d'))))
    
            % Base polygon and z-range
            Ra = 1;                                     % Polygon apothem
            zmin = min(V(:,3));                         % Minimum z-point
            zmax = max(V(:,3));                         % Maximum z-point
            Vxy = V(V(:,3)==zmin,[1,2]);                % Vertices on the xy plane
            thVxy = atan2(Vxy(:,2),Vxy(:,1));           % Angle of vertices
            thVxy(thVxy<0) = thVxy(thVxy<0)+2*pi;
            [thVxy,Ixy] = sort(thVxy);                  % Sorting
            Vxy = Vxy(Ixy,:);
            if abs(thmax-2*pi)<=2*eps                
                nV = size(Vxy,1)/4;                     % Number of sides per quadrant       
                Vxy = Vxy([end,1:end,1],:);             % To allow for interpolation past 2*pi
                thVxy = [(thVxy(end)-2*pi);thVxy;(thVxy(1)+2*pi)];
            else
                nV = (size(Vxy,1)-2)/2;                       
            end
            dthVxy = pi/(2*nV);                         % Step in theta of the polygon
            
            % Factors in thvmax steps to preserve homogeneous sampling density
            thvfactor = thmax/(pi/2);
    
            % Sampling
            [R,th] = deal(nan(max(1,4*steps^2),1));                 % Preallocating
            Rvth = 0:(pi/(2*steps)):(pi/2);                         % Linear sampling on 2D rotation angle
            Rv = tan(Rvth/2);                                       % Irregular sampling on 2D Rodrigues vector
            c = 0;                                                  % Counter
            for j = 1:length(Rv)
                % Steps in th
                if j < length(Rv)
                    stepsth = steps;
                else                                                % To preferentially sample corners
                    stepsth = round(nV*round(steps/nV));            % Factor
                end
                steps_thvj = ceil(Rv(j)*stepsth);                   % Number of steps proportial to the circumference
                thv = 0:thmax/ceil(thvfactor*steps_thvj):thmax;
                if length(thv) <= 1                                 % At least point
                    thv = 0;
                end

                % Maximum distance of the polygon as a function of th
                thetaa = mod(thv+dthVxy/2,dthVxy)-dthVxy/2;         % Angle between thv and the closest apothem
                Rmaxth = Ra./cos(thetaa);                           % Maximum distance in thv direction
                % Steps in th
                for k = 1:length(thv)
                    c = c+1;
                    R(c) = Rv(j)*Rmaxth(k);                         % r weighted by Rmaxth
                    th(c) = thv(k);                                 % Symmetry range
                end
            end
            R(c+1:end) = [];                                        % Trimming
            th(c+1:end)= [];
    
            % Extruding in z
            zmaxth = 2*atan(zmax);                                  % Maximum rotation in z
            zminth = 2*atan(zmin);
            zsteps = round(zmax*steps);
            dzth = zmaxth/zsteps;                                   % Step in z
            if steps >= 1
                dzth = min(dzth,zmaxth);
            end
            zvth = (zminth:dzth:zmaxth)';                           % Linear sampling on rotation rz
            zv = tan(zvth/2);                                       % Irregular sampling on rz Rodrigues vector                                 
            
            R  = repmat( R,length(zv),1);                           % Repeating elements
            th = repmat(th,length(zv),1);
            rz = repelem(zv,c);
    
            % Change to cartesian coordinates
            rv = [R.*cos(th),R.*sin(th),rz];                        % Rodrigues vectors
    
        
        %%%%% D2h-[D2h,C2h,Ci] - square-prism grid
        elseif (ge == 0) && ...
               ((strcmp(pg2,'D2h') && (strcmp(pg1,'D2h') || strcmp(pg1,'C2h') || strcmp(pg1,'C2h*') || strcmp(pg1,'Ci'))))   
        
            % Minimum & maxima
            xminth = 2*atan(min(V(:,1)));                         
            xmaxth = 2*atan(max(V(:,1)));
            yminth = 2*atan(min(V(:,2)));
            ymaxth = 2*atan(max(V(:,2)));
            zminth = 2*atan(min(V(:,3)));
            zmaxth = 2*atan(max(V(:,3)));
        
            % Grid in each dimension
            dth = pi/(2*steps);                                     % Step in z
            xvth = (xminth:dth:xmaxth)';                            % Linear sampling on rotation angle
            xv = tan(xvth/2);                                       % Irregular sampling on rz Rodrigues vector
            yvth = (yminth:dth:ymaxth)';                            % Linear sampling on rotation angle
            yv = tan(yvth/2); 
            zvth = (zminth:dth:zmaxth)';
            zv = tan(zvth/2); 
            [rx,ry,rz] = meshgrid(xv,yv,zv);                        % Total grid        
            
            % Rodrigues vectors
            rv = nan(numel(rx),3);                                  
            rv(:,1) = rx(:);
            rv(:,2) = ry(:);
            rv(:,3) = rz(:);
    
        %%%%% Not coded
        else
            error(['The fundamental zone ',pg1,'-',pg2,' grid is not included.'])
        end


    %%%%%%%%%%%%%%%%% UNIQUE
    elseif(strcmp(steps,'unique'))
        %%%%% Oh-Oh
        if (ge == 1) && strcmp(pg1,'Oh') && strcmp(pg2,'Oh') 
            rv1max = sqrt(2)-1;                 % Largest rodrigues vector component
            c = 0;                              % Counter
            rv = nan(21,3);                     % Preallocating
            
            % Rodrigues vectors
    %         c=c+1; % O
    %         rv(c,:) = [0,0,0];          
            c=c+1; % A
            rv(c,:) = [rv1max,0,0];     
            c=c+1; % B
            rv(c,:) = [rv1max,rv1max,0];
            c=c+1; % E
            rv(c,:) = [1/3,1/3,1/3];
            c=c+1; % C
            rv(c,:) = [rv1max,rv1max,1-2*rv1max];
            c=c+1; % D
            rv(c,:) = [rv1max,1/2-rv1max/2,1/2-rv1max/2];
            c=c+1; % OA
            rv(c,:) = 1/2*[rv1max,0,0];
            c=c+1; % OB
            rv(c,:) = 1/2*[rv1max,rv1max,0];
            c=c+1; % OE
            rv(c,:) = 1/2*[1/3,1/3,1/3];
            c=c+1; % AB
            rv(c,:) = [rv1max,rv1max/2,0];
            c=c+1; % AD
            rv(c,:) = 1/2*([rv1max,0,0] + [rv1max,1/2-rv1max/2,1/2-rv1max/2]);
            c=c+1; % BC
            rv(c,:) = [rv1max,rv1max,1/2*(1-2*rv1max)];
            c=c+1; % CD
            rv(c,:) = 1/2*([rv1max,rv1max,1-2*rv1max] + [rv1max,1/2-rv1max/2,1/2-rv1max/2]);
            c=c+1; % CE
            rv(c,:) = 1/2*([rv1max,rv1max,1-2*rv1max] + [1/3,1/3,1/3]);
            c=c+1; % DE
            rv(c,:) = 1/2*([rv1max,1/2-rv1max/2,1/2-rv1max/2] + [1/3,1/3,1/3]);
            c=c+1; % *AC
            rv(c,:) = 1/2*([rv1max,0,0] + [rv1max,rv1max,1-2*rv1max]);
            c=c+1; % OAB
            rv(c,:) = 1/3*([rv1max,0,0] + [rv1max,rv1max,0]);
            c=c+1; % OBCE
            rv(c,:) = 1/4*([rv1max,rv1max,0] + [rv1max,rv1max,1-2*rv1max] + [1/3,1/3,1/3]);
            c=c+1; % OADE
            rv(c,:) = 1/4*([rv1max,0,0] + [rv1max,1/2-rv1max/2,1/2-rv1max/2] + [1/3,1/3,1/3]);
            c=c+1; % ABCD
            rv(c,:) = 1/4*([rv1max,0,0] + [rv1max,rv1max,0] + [rv1max,rv1max,1-2*rv1max] + [rv1max,1/2-rv1max/2,1/2-rv1max/2])+[0,0,0.1];
            c=c+1; % CDE
            rv(c,:) = 1/3*([rv1max,rv1max,1-2*rv1max] + [rv1max,1/2-rv1max/2,1/2-rv1max/2] + [1/3,1/3,1/3]);
            c=c+1; % interior
            rv(c,:) = 1/4*([rv1max,0,0] + [rv1max,rv1max,0] + [1/3,1/3,1/3]);

        %%%%% D2h-[D4h,D6h,D8h,D3d] - Polygon extruded in z
        elseif (ge == 0) && ...
               ((strcmp(pg1,'D2h') && (strcmp(pg2,'D8h') || strcmp(pg2,'D6h') || strcmp(pg2,'D4h') || strcmp(pg2,'D3d')))) 
        % TODO
    
        %%%%% D2h-[D2h,C2h,Ci] - square-prism grid
        elseif (ge == 0) && ...
               ((strcmp(pg2,'D2h') && (strcmp(pg1,'D2h') || strcmp(pg1,'C2h') || strcmp(pg1,'C2h*') || strcmp(pg1,'Ci'))))
        % TODO
    
        %%%%% Not coded
        else
            error(['The fundamental zone ',pg1,'-',pg2,' unique points are not included.'])
            
        end
    end

    %%%%%%%%%%%%%%%%% REGION
    if reg == 1
        [region,w] = regionSO3(pg1,pg2,rv,ge);
    else
        [region,w] = deal(nan);
    end
    
   
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plots

%%% To visualize the FZ:
% figure
% for i = 1:length(F)                   % Plot FZ
%     plygn = fill3(V(F{i},1),V(F{i},2),V(F{i},3),'r','EdgeColor','r','LineWidth',3);
%     hold on                
%     plygn.FaceAlpha = 0.1;
% end
% scatter3(rv(:,1),rv(:,2),rv(:,3),10,'b','filled')    % Grid
% axis equal
% xlabel('r_x')
% ylabel('r_y')
% zlabel('r_z')
% set(gcf,'color','w')
% grid on
%

