%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the region within an SO(3) fundamental zone a point 
% belongs to.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Inputs:
% pg1,pg2           Point groups of the two elements 
% rv                Rodrigues vector of each point
% ge                Grain exchange symmetry? 0 -> no, 1 -> yes
%
%%% Outputs: 
%   region          Region within the FZ:
%   w               Weights (i.e. number of symmetric points within the FZ)
%                   To be implemented for multiple FZ geometries. 
%

function [region,w] = regionSO3(pg1,pg2,rv,ge)

    pgsorted = sort(string({pg1,pg2}));         % Sort point groups
    pg1 = pgsorted(1);
    pg2 = pgsorted(2);

    r1 = rv(:,1);                               % Rodrigues vector components
    r2 = rv(:,2);
    r3 = rv(:,3);

    %%%%% Oh-Oh
    if (ge == 1) && strcmp(pg1,'Oh') && strcmp(pg2,'Oh') 



        t = 20*eps;                      % Tolerance
        rmax = sqrt(2)-1;

        region = cell(size(rv,1),1);    % Preallocation
        w = nan(size(rv,1),1);

        % Points
        id = abs(r1)<=t & abs(r2)<=t & abs(r3)<=t & isnan(w);
            region(id) = {'O'};
            w(id) = 1;
        id = abs(r1-rmax)<=t & abs(r2)<=t & abs(r3)<=t & isnan(w);
            region(id) = {'A'};
            w(id) = 6;
        id = abs(r1-rmax)<=t & abs(r2-rmax)<=t & abs(r3)<=t & isnan(w);
            region(id) = {'B'};
            w(id) = 12;
        id = abs(r1-rmax)<=t & abs(r2-rmax)<=t & abs(r3-1+2*rmax)<=t & isnan(w);
            region(id) = {'C'};
            w(id) = 24;
        id = abs(r1-rmax)<=t & abs(r2-1/2+rmax/2)<=t & abs(r3-1/2+rmax/2)<=t & isnan(w);
            region(id) = {'D'};       % *every point D is also a point B
            w(id) = 24;
        id = abs(r1-1/3)<=t & abs(r2-1/3)<=t & abs(r3-1/3)<=t & isnan(w);
            region(id) = {'E'};
            w(id) = 8;
    
        % Lines
        id = abs(r2)<=t & abs(r3)<=t & isnan(w);
            region(id) = {'OA'};
            w(id) = 6;
        id = abs(r1-r2)<=t & abs(r3)<=t & isnan(w);
            region(id) = {'OB'};
            w(id) = 12;
        id = abs(r1-r3)<=t & abs(r2-r3)<=t & isnan(w);
            region(id) = {'OE'};
            w(id) = 8;
        id = abs(r1-rmax)<=t & abs(r3)<=t & isnan(w);
            region(id) = {'AB'};
            w(id) = 24;
        id = abs(r1-rmax)<=t & abs(r2-r3)<=t & isnan(w);
            region(id) = {'AD'};      % *every point AD is also a point AB
            w(id) = 24;
        id = abs(r1-rmax)<=t & abs(r2-rmax)<=t & isnan(w);
            region(id) = {'BC'};
            w(id) = 24;
        id = abs(r1-rmax)<=t & abs(r2+r3-1+rmax)<=t & isnan(w);
            region(id) = {'CD'};      % *every point CD is also a point BC
            w(id) = 48;
        id = abs(r1+r3/2-1/2)<=t & abs(r2+r3/2-1/2)<=t & isnan(w);
            region(id) = {'CE'};
            w(id) = 24;
        id = abs(r1+2*r3-1)<=t & abs(r2-r3)<=t & isnan(w);
            region(id) = {'DE'};
            w(id) = 24;
        id = abs(r1-rmax)<=t & abs(r2+rmax*r3/(2*rmax-1))<=t & isnan(w);
            region(id) = {'*AC'};
            w(id) = 48;
    
        % Surfaces
        id = abs(r3)<=t & isnan(w);
            region(id) = {'OAB'};
            w(id) = 24;
        id = abs(r1-rmax)<=t & isnan(w);
            region(id) = {'ABCD'};
            w(id) = 48;
        id = abs(r1+r2+r3-1)<=t & isnan(w);
            region(id) = {'CDE'};
            w(id) = 48;
        id = abs(r2-r1)<=t & isnan(w);
            region(id) = {'OBCE'};
            w(id) = 24;
        id = abs(r3-r2)<=t & isnan(w);
            region(id) = {'OADE'};
            w(id) = 24;

        % Interior
        id = isnan(w);
            region(id) = {'interior'};
            w(id) = 24;

    %%%%% D2h-[D4h,D6h,D8h,D3d,D2h,C2h,Ci] - Polygon extruded in z 
    % TODO: weights
    elseif (ge == 0) && (((strcmp(pg1,'D2h') && (strcmp(pg2,'D8h') || strcmp(pg2,'D6h') || strcmp(pg2,'D4h') || strcmp(pg2,'D3d')))) || ... 
                          (strcmp(pg2,'D2h') && (strcmp(pg1,'D2h') || strcmp(pg1,'C2h') || strcmp(pg1,'C2h*') || strcmp(pg1,'Ci'))))
    
        % Fundamental zone in Rodrigues-Frank space
        [V,F,thmax] = fzSO3(pg1,pg2);               % Vertices and faces

        % Base polygon and z-range
        zminmax = [min(V(:,3)),max(V(:,3))];        % Minimum and maximum z-coordinates
        Vxy = V(V(:,3)==zminmax(1),[1,2]);          % Vertices on the xy plane
        thVxy = atan2(Vxy(:,2),Vxy(:,1));           % Angle of vertices
        thVxy(thVxy<0) = thVxy(thVxy<0)+2*pi;
        [~,Ixy] = sort(thVxy);                  % Sorting
        Vxy = Vxy(Ixy,:);
        alphabet = 'ABCDEFGHIJKLMNPQRSTUVWXYZ';     % Assigning point letters (O first)
        points = alphabet(1:size(Vxy,1))';
        points(sum(Vxy==0,2)==2) = 'O';             % If point O is a vertex
        if strcmp(pg1,'C2h*')                       % Adjusting points for C2h*
            points = points([3,4,1,2]);
        end

        % FZ size within the whole orientation space
        nth = round(2*pi/thmax);                    % Times the FZ is repeated in theta
        nz = round(pi/diff(atan(zminmax)*2));       % Times the FZ is repeated in z

        % Assigning
        t = 10*eps;                                 % Tolerance
        region = cell(size(rv,1),1);                % Preallocation
        w = nan(size(rv,1),1);                      

        % Points
        id = abs(r1)<=t & abs(r2)<=t & abs(r3)<=t & isnan(w);   % O
            region(id) = {'O'};
            w(id) = round(nz/2);
        for i = 1:length(points)                    % Lower surface
            if zminmax(1) < 0 
                lbl1 = 'm';  % minus
            else
                lbl1 = '';   % z = 0
            end
            id = abs(r1-Vxy(i,1))<=t & abs(r2-Vxy(i,2))<=t & abs(r3-zminmax(1))<=t & isnan(w);
            region(id) = {[points(i),lbl1]};
            w(id) = nth*nz/2;   %TODO
        end
        for i = 1:length(points)                    % Upper surface
            lbl2 = 'p';      % plus
            id = abs(r1-Vxy(i,1))<=t & abs(r2-Vxy(i,2))<=t & abs(r3-zminmax(2))<=t & isnan(w);
            region(id) = {[points(i),lbl2]};
            w(id) = nth*nz/2;   %TODO
        end

        % Edges
        Vxye = [Vxy(2:end,:);Vxy(1,:)];
        pointse = [points(2:end);points(1)];
        for i = 1:length(points)                    % Vertical lines
            id = abs(r1-Vxy(i,1))<=t & abs(r2-Vxy(i,2))<=t & isnan(w);
            region(id) = {[points(i),lbl1,points(i),lbl2]};
            w(id) = nth*nz/2;   %TODO
        end
        for i = 1:length(points)                    % Lower surface
            id = dist2point([r1,r2],Vxy(i,:),Vxye(i,:))<=t & abs(r3-zminmax(1))<=t & isnan(w);
            region(id) = {[points(i),lbl1,pointse(i),lbl1]};
            w(id) = nth*nz/2;   %TODO
        end
        for i = 1:length(points)                    % Upper surface
            id = dist2point([r1,r2],Vxy(i,:),Vxye(i,:))<=t & abs(r3-zminmax(2))<=t & isnan(w);
            region(id) = {[points(i),lbl2,pointse(i),lbl2]};
            w(id) = nth*nz/2;   %TODO
        end

        % Faces
        for i = 1:length(points)                    % Vertical faces
            id = dist2point([r1,r2],Vxy(i,:),Vxye(i,:))<=t & isnan(w);
            region(id) = {[points(i),lbl1,pointse(i),lbl1,pointse(i),lbl2,points(i),lbl2]};
            if nth == 2 && i == length(points)              % Long vertical face
                w(id) = nz;     %TODO
            else                                            % All other vertical faces (it may actually vary depending on face direction)
                w(id) = nth*nz; %TODO
            end
        end
        id = abs(r3-zminmax(1))<=t & isnan(w);      % Lower surface
            reg = [points,repmat(lbl1,length(points),1)]';
            region(id) = {reg(:)'};
            w(id) = nth*nz/2;
        id = abs(r3-zminmax(2))<=t & isnan(w);      % Upper surface
            reg = [points,repmat(lbl2,length(points),1)]';
            region(id) = {reg(:)'};
            w(id) = nth*nz/2;

        % Interior
        id = isnan(w);
            region(id) = {'interior'};
            w(id) = nth*nz;
            
    %%%%% Not coded
    else
        error(['The fundamental zone ',pg1,'-',pg2,' is not included.'])

    end

end

% Function: distance from points to a line define by vertices v1 & v2
function d = dist2point(points,v1,v2)

    % Ax + By + C = 0.
    A = v2(2) - v1(2);
    B = v1(1) - v2(1);
    C = v2(1)*v1(2) - v1(1)*v2(2);
      
    % Calculate distances as d = |Ax+By+C|/sqrt(A^2+B^2)
    d = abs(A*points(:, 1)+B*points(:, 2) + C) / sqrt(A^2+B^2);
end



%%%%% Plots

%%% To visualize the FZ:
% figure
% for i = 1:length(F)                             % Plot FZ
%     plygn = fill3(V(F{i},1),V(F{i},2),V(F{i},3),'r','EdgeColor','r','LineWidth',3,'HandleVisibility','off');
%     hold on                
%     plygn.FaceAlpha = 0.1;
% end
% ureg = unique(region);                          % Assigning colors
% clrs = distinguishable_colors(length(ureg));
% clr = nan(length(ureg),1);
% for i = 1:length(ureg)
%     id = cellfun(@(x) strcmp(x,ureg{i}),region);
%     scatter3(rv(id,1),rv(id,2),rv(id,3),10,clrs(i,:),'filled','DisplayName',ureg{i})    % Grid
% end
% axis equal
% xlabel('r_x')
% ylabel('r_y')
% zlabel('r_z')
% legend(ureg,'Location','EastOutside','FontSize',min(10,200/length(ureg)))
% set(gcf,'color','w')
% grid on
%