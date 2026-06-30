%%% Example - Damage distributions for uniaxial load of an fcc metal

% Inputs                         
    % Material - 304 stainless steel
s11 = 1.22e-2;                      % Compliance components
s12 = -5.02e-3;
s44 = 9.12e-3;
fI = 0.5;                           % Volume fraction of grain I
type_load = 'uniaxial';             % Uniaxial stress
stress = 1;                         % Stress magnitude
steps_dis = 4;                      % Steps in disorientation FZ
steps_bp = 8;                       % Steps in boundary plane FZ
steps_load = 8;                     % Steps in stress FZ
outputs = {'mS','strainenergy'};    % Outputs from the elastic bicrystal model


%%%    ---------    Analysis    ---------    %%%

% Grids
dis = grid_disorientation(steps_dis);                   % Disorientation 
gb = grid_boundaryplane(dis,steps_bp);                  % Boundary plane
tests = grid_stress(gb,type_load,stress,steps_load);    % Stress

% Incompatibility stress
Ntests = size(tests.id,1);                              % Number of conditions
mSmax = nan(Ntests,1);                                  % Preallocating
for i = 1:Ntests
    % Elastic bicrystal incompatibility stress model
    [ssI,ssII,stressCB] = bicrystal([s11,s12,s44],fI,tests.stressR(i,:)',tests.euleranglesI(i,:),tests.euleranglesII(i,:),tests.GBangles(i,:),outputs);

    % Calculating damage metrics
    mSmax(i) = max(abs([ssI.mS*ssI.stressmag;ssII.mS*ssII.stressmag]))/stress;  % Schmid factors
end

%%%    ---------    Weighted distributions & plots    ---------    %%%

% Statistics
mSmax_stats = symmetryresults(mSmax,tests.id,0,'stats',tests.w);
disp('mSmax weighted statistics:')
disp(mSmax_stats)

% Overall distribution
Nbins = 20;                                                     % Number of bins
Wi = tests.w(:,1).*tests.w(:,2).*tests.w(:,3);                  % Total weight for each test point
Wtotal = sum(Wi);                                               % Total number of points accounting for the symmetries
edges = min(mSmax):(max(mSmax)-min(mSmax))/Nbins:max(mSmax);    % Bin edges
binsize = edges(2)-edges(1);                                    % Bin size
bincenters = (edges(1:end-1)+edges(2:end))/2;                   % Bin centers
bin = arrayfun(@(x) find(edges<=x,1,'last'),mSmax);             % Bin each condition belongs to
mSmaxhist = nan(1,length(edges)-1);                             % Preallocating
for i = 1:(length(edges)-1)
    mSmaxhist(i) = sum(Wi(bin == i));                           % Sum points per bin
end

figure
plot(bincenters,mSmaxhist/(Wtotal*binsize),'LineWidth',1.5)
xlabel('$m_{GB,\textrm{max}}$','Interpreter','Latex')
ylabel('Probability density')
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',12)
box on

% Conglomerating data by disorientation type
fzu = unique(tests.fzdichromatic);                  % unique dichromatic FZ regions
IDdis = tests.id;
for i = 1:length(fzu)
    idu = strcmp(tests.fzdichromatic,fzu{i});       % conditions for each region
    IDdis(idu,1) = i;                               % grouping
end
[mSmax_s,rows] = symmetryresults(mSmax,IDdis,1,'stats',tests.w);          % Conglomerating results
tests_s = sh(tests,rows); 

% Conglomerating data by disorientation angle theta
dtheta = 10;                                % Step size in theta
thetaedges = 0:dtheta:70;                   % theta bins
theta = tests.th*180/pi;                    % theta angle
mSmax_theta = binquart(mSmax,theta,thetaedges,tests.w); % Binning function

figure
boxplot(mSmax_theta',thetaedges(2:end),'Whisker',inf,'Width',0.8)
xticks((1:length(thetaedges))-0.5)
xticklabels(thetaedges)
xlabel('$\theta \, [^{\circ}]$','Interpreter','Latex')
ylabel('$m_{GB,\textrm{max}}$','Interpreter','Latex')
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',12)
box on

% Conglomerating data by angle alpha (between disorientation axis & uniaxial load)
dalpha = 10;                                % Step size in alpha
alphaedges = 0:dalpha:90;                   % alpha bins
alphaedges(end) = 90.00000001;              % So that last bin includes 90 deg
alpha = nan(size(tests.id,1),1);            % Preallocating
for i = 1:size(tests.id,1)                  % Calculating alpha
    alpha(i) = acosd(dot(tests.l(i,:),tests.nloadaxis(i,:)));
    if alpha(i) > 90
        alpha(i) = 180 - alpha(i);
    end
end
mSmax_alpha = binquart(mSmax,alpha,alphaedges,tests.w); % Binning function

figure
boxplot(mSmax_alpha',alphaedges(2:end),'Whisker',inf,'Width',0.8)
xticks((1:length(alphaedges))-0.5)
xticklabels(alphaedges)
xlabel('$\alpha \, [^{\circ}]$','Interpreter','Latex')
ylabel('$m_{GB,\textrm{max}}$','Interpreter','Latex')
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',12)
box on

% Conglomerating data by angle phi (between boundary plane normal & uniaxial load)
dphi = 10;                                  % Step size in phi
phiedges = 0:dphi:90;                       % phi bins
phiedges(end) = 90.00000001;                % So that last bin includes 90 deg
phi = nan(size(tests.id,1),1);              % Preallocating
for i = 1:size(tests.id,1)                  % Calculating phi
    phi(i) = acosd(dot(tests.nboundary(i,:),tests.nloadaxis(i,:)));
    if phi(i) > 90
        phi(i) = 180 - phi(i);
    end
end
mSmax_phi = binquart(mSmax,phi,phiedges,tests.w);

figure
boxplot(mSmax_phi',phiedges(2:end),'Whisker',inf,'Width',0.8)
xticks((1:length(phiedges))-0.5)
xticklabels(phiedges)
xlabel('$\alpha \, [^{\circ}]$','Interpreter','Latex')
ylabel('$m_{GB,\textrm{max}}$','Interpreter','Latex')
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',12)
box on

