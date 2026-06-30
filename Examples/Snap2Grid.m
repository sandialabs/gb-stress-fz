%%% Example - Mapping stressed grain boundary configurations to the FZ
% Example runs with 500 randomly generated stressed bicrystals, and 
% mapping is performed snapping the configurations to their closest points
% within a predefined grid.

% Inputs                         
    % Material - 304 stainless steel
s11 = 1.22e-2;                      % Compliance components
s12 = -5.02e-3;
s44 = 9.12e-3;
fI = 0.5;                           % Volume fraction of grain I
outputs = 'mS';                     % Outputs from the elastic bicrystal model
rng(12345);                         % Setting an RNG seed 

                                    % Randomly generated stressed bicrystals
Ntests = 1000;                                          % Number of configurations
euleranglesI  = randv.randEuler(Ntests);                % Euler angles grain I
euleranglesII = randv.randEuler(Ntests);                % Euler angles grain II
GBangles = randv.randVector(Ntests,'vectorangles');     % Grain boundary angles
stress = randv.randStress(1,'uniaxial');                % Stress with norm of 1

% Predefined grain boundary grid
% Option 'includeedges' to make sure points repeated points at the edges of
% the disorientation  and boundary plane grids are kept
dis_grid = grid_disorientation(4);
gb_grid = grid_boundaryplane(dis_grid,8,'symmetry','includeedges'); 

% Map to FZ using 'snap2grid' option:
% This applies for functions 'toFZ_boundaryplane' too.
% For easier data visualization try: struct2table(tests)
tests = toFZ_stress(euleranglesI,euleranglesII,GBangles,stress,'snap2grid',gb_grid);

disp('Angles between the original and the grid points:')
disp(['- Minimum: ',num2str(min(tests.th_gbgrid)),' °'])
disp(['- Maximum: ',num2str(max(tests.th_gbgrid)),' °'])
disp(['- Mean: ',num2str(mean(tests.th_gbgrid)),' ± ',num2str(std(tests.th_gbgrid)),' °'])


% Incompatibility stress - used to proove that the original and mapped
% configurations are indeed identical
[rssmax,rssmaxFZ] = deal(nan(Ntests,1));                  % Preallocating
for i = 1:Ntests
    % Original
    [ssI,ssII,~] = bicrystal([s11,s12,s44],fI,stress',euleranglesI(i,:),euleranglesII(i,:),GBangles(i,:),outputs);
    rssmax(i) = max(abs([ssI.mS*ssI.stressmag;ssII.mS*ssII.stressmag]));    % Highest resolved shear stress

    % FZ equivalent
    [ssI,ssII,~] = bicrystal([s11,s12,s44],fI,tests.stressR(i,:)',tests.euleranglesI(i,:),tests.euleranglesII(i,:),tests.GBangles(i,:),outputs);
    rssmaxFZ(i) = max(abs([ssI.mS*ssI.stressmag;ssII.mS*ssII.stressmag]));  % Highest resolved shear stress
end

%%% Plots
fig = figure;                       % Comparing original to FZ 
scatter(rssmax,rssmaxFZ,10,'filled')
hold on
plot([0,1],[0,1],'--','Color','r')
xlabel('$RSS_{\textrm{max}}$','Interpreter','latex')
ylabel('$RSS_{\textrm{max,snapped}}$','Interpreter','latex')
set(gca,'FontName','Arial','FontSize',14)
axis equal
xlim([0,0.7])
ylim([0,0.7])
set(gcf,'color','w');
grid off
box on