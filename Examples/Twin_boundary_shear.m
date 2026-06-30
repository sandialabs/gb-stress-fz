%%% Example - Pure shear on an fcc coherent twin boundary

% Inputs                         
    % Material - 304 stainless steel
s11 = 1.22e-2;                      % Compliance components
s12 = -5.02e-3;
s44 = 9.12e-3;
fI = 0.5;                           % Volume fraction of grain I
type_dis = 'twin';                  % Sigma3 twin misorientation
type_load = 'shear';                % Shear stress
stress = 1;                         % Stress magnitude
steps_bp = 0;                       % Steps in boundary plane FZ
steps_load = 90;                    % Steps in stress FZ
outputs = {'mS','strainenergy'};    % Outputs from the elastic bicrystal model


%%%    ---------    Analysis    ---------    %%%

% Grids
dis = grid_disorientation(type_dis);                    % Disorientation 
gb = grid_boundaryplane(dis,steps_bp);                  % Boundary plane
tests = grid_stress(gb,type_load,stress,steps_load);    % Stress

% Incompatibility stress
Ntests = size(tests.id,1);                              % Number of conditions
[mSmax,Rgrain,strainenergy,] = deal(nan(Ntests,1));     % Preallocating
for i = 1:Ntests
    % Elastic bicrystal incompatibility stress model
    [ssI,ssII,stressCB] = bicrystal([s11,s12,s44],fI,tests.stressR(i,:)',tests.euleranglesI(i,:),tests.euleranglesII(i,:),tests.GBangles(i,:),outputs);

    % Calculating damage metrics
    mSmax(i) = max(abs([ssI.mS*ssI.stressmag;ssII.mS*ssII.stressmag]))/stress;  % Schmid factors
    mSIsorted = sort(abs(ssI.mS*ssI.stressmag),'descend')/stress;               % Sorted Schmid factors for each grain
    mSIIsorted = sort(abs(ssII.mS*ssII.stressmag),'descend')/stress;  
    RB = [mSIsorted(1)/mSIsorted(2),mSIIsorted(1)/mSIIsorted(2)];               % R ratio for each grain
    if abs(mSIsorted(1)-mSIIsorted(1))<10*eps
        Rgrain(i) = max(RB);                                                    % If both grains active, select highest R ...
    else
        [~,IRmax] = max([mSIsorted(1),mSIIsorted(1)]);  
        Rgrain(i) = RB(IRmax);                                                  % ... otherwise R for the grain with the highest mS
    end
    strainenergy(i) = (ssI.strainenergy+ssII.strainenergy)/(2*s11*stress^2);    % Normalized mean strain energy density
end

%%%    ---------    Plots    ---------    %%%

% m_{gb,max} highest Schmid factor
fig = plot_stress(tests,'all',mSmax);
colorbar
fig{1}.Children(1).Label.String = '$m_{GB,\textrm{max}}$';
fig{1}.Children(1).Label.Interpreter = 'latex';
    
% R_{GB} ratio between the highest and second highest Schmid factor
fig = plot_stress(tests,'all',Rgrain);
colorbar
fig{1}.Children(1).Label.String = '$R_{GB}$';
fig{1}.Children(1).Label.Interpreter = 'latex';

% u_{GB} normalized mean strain energy density
fig = plot_stress(tests,'all',strainenergy);
colorbar
fig{1}.Children(1).Label.String = '$u_{GB}$';
fig{1}.Children(1).Label.Interpreter = 'latex';

% m_{gb,max} vs R_{GB} pairs
fig = figure;
scatter(mSmax,Rgrain,10,'filled')
xlabel('$m_{GB,\textrm{max}}$','Interpreter','latex')
ylabel('$R_{GB}$','Interpreter','latex')
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',12)
box on

