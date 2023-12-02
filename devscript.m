mrstModule add baia-mare
mrstModule add test-suite
mrstModule add ad-core ad-props compositional
mrstModule add geothermal
mrstModule add upr
mrstModule add optimization
mrstModule add mrst-gui

mrstVerbose on

loadData = @(name) load(fullfile(mrstPath('baia-mare'), name));

%% Set up test case
setup = TestCase('baia_mare_geothermal_2', 'readGridFromDisk', true);
setup.plot();

%% Set up and simulate problem
problem = setup.getPackedSimulationProblem();
simulatePackedProblem(problem, 'restartStep', nan);

%% Get results
[~, states, reports] = getPackedSimulatorOutput(problem);
setup.plot(states);

%% Set up optimization problem (TBC)
pv = setup.model.operators.pv;
% set up 'matrix' for parameter options for easier editing. The specific
% limits set for the various parameters influences the tuning/optimization % procedure to a large extent

G = setup.model.G;
pcell = G.cells.tag.crus_cell.*1 + G.cells.tag.sed1_cell.*2 + G.cells.tag.sed2_cell.*3;
pface = max(pcell(setup.model.operators.N), [], 2);

config = {...
     %name           include    scaling              boxlims   relativeLimits  
    'porevolume',       1,     'linear',    [.01*pv, 1.5*pv],              [], pcell
    'transmissibility', 1,        'log',                  [],     [1e-2  1e2], pface
};

parameters = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    parameters = addParameter(parameters, setup, ...
        'name',    config{k,1}, 'scaling', config{k,3}, ...
        'boxLims', config{k,4}, 'relativeLimits',config{k,5});
end

%% Model calibration
statesRef = states;
pvec = getScaledParameterVector(setup, parameters);
mismatchFn = @(model, states, schedule, states_ref, compDer, tstep, state) ...
    temperatureMismatch(model, states, schedule, states_ref);
objh2 = @(p) evaluateMatch(p, mismatchFn, setup, parameters, statesRef);
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 30 iterations
[v2, p_opt2, history2] = unitBoxBFGS(pvec, objh2, 'maxIt', 20);
               
%%
FF = scatteredInterpolant;
FF.Points = setup.model.G.cells.centroids; 
FF.ExtrapolationMethod='none';
FF.Values = states{end}.T;

%%
DemWelTVal = loadData('DemWell1-6.mat');
data1 = loadData('Wells1.mat'); Wells = data1.Wells;
data2 = loadData('Wells2.mat'); DemWel = data2.DemWel;

figure(f),  hold on,

plot(Wells.profdstA05,Wells.Elevation,'^','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10)
text(Wells.profdstA05,Wells.Elevation,Wells.Var3)

plot(DemWel.profdstA05,DemWel.Elevation,'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10)
text(DemWel.profdstA05,DemWel.Elevation,DemWel.Site), ylim([-9e3 2e3])

%%
n = 100;

% z = DemWelTVal.well2(:,2)-DemWel.Elevation(4);

z = -[0 1000 2000 2500]';
x = repmat(DemWel.profdstA05(5), numel(z), 1);  

x = [x, z];
cells = findEnclosingCell(setup.model.G, x);

close all
plotGrid(setup.model.G, cells);
plotGrid(setup.model.G, 'FaceColor', 'none');
setup.setAxisProperties(gca);


%%

zz = linspace(-1500,600,20)';

figure(40),clf

subplot(2,4,1), plot(DemWelTVal.well2(:,1),DemWelTVal.well2(:,2)-DemWel.Elevation(2),'b'),set(gca,'Ydir','reverse'), title(DemWel.Site(2)),%ylim([-1000 1000]),xlim([10 30])
tt = FF(DemWel.profdstA05(2) + zz*0, zz); 
hold on, plot(tt,-zz,'r'), legend('Data', 'Model')

subplot(2,4,2), plot(DemWelTVal.well3(:,1),DemWelTVal.well3(:,2)-DemWel.Elevation(3),'b'),set(gca,'Ydir','reverse'), title(DemWel.Site(3)),ylim([-1000 1000]),xlim([10 30])
tt = FF(DemWel.profdstA05(3) + zz*0, zz); 
hold on, plot(tt,-zz,'r'), legend('Data', 'Model')

subplot(2,4,3), plot([Wells.zo(2),Wells.z1000(2),Wells.z2000(2),Wells.z2500(2)],...
    [0 1000 2000 2500],'b'),set(gca,'Ydir','reverse'), title(Wells.Var3(2)),%ylim([-1000 1000]),xlim([10 150])
tt = FF(Wells.profdstA05(2) + zz*0, zz); 
hold on, plot(tt,-zz,'r'), legend('Data', 'Model')

subplot(2,4,4), plot([Wells.zo(1),Wells.z1000(1),Wells.z2000(1),Wells.z2500(1)],...
    [0 1000 2000 2500],'b'),set(gca,'Ydir','reverse'), title(Wells.Var3(1)),%ylim([-1000 1000]),xlim([10 150])
hold on, plot(DemWelTVal.well6(:,1),DemWelTVal.well6(:,2)-DemWel.Elevation(6),'b'),
tt = FF(Wells.profdstA05(1) + zz*0, zz); 
hold on, plot(tt,-zz,'r'), legend('Data', 'Model')

set(gca,'Ydir','reverse'), title(DemWel.Site(6)),%ylim([-1000 1000]),xlim([10 150])
subplot(2,4,5), plot([Wells.zo(6),Wells.z1000(6),Wells.z2000(6),Wells.z2500(6)],...
    [0 1000 2000 2500],'b'),set(gca,'Ydir','reverse'), title(Wells.Var3(6)),%ylim([-1000 1000]),xlim([10 150])
tt = FF(Wells.profdstA05(6) + zz*0, zz); 
hold on, plot(tt,-zz,'r'), legend('Data', 'Model')

subplot(2,4,6), plot([Wells.zo(3),Wells.z1000(3),Wells.z2000(3),Wells.z2500(3)],...
    [0 1000 2000 2500],'b'),set(gca,'Ydir','reverse'), title(Wells.Var3(3)),%ylim([-1000 1000]),xlim([10 150])
tt = FF(Wells.profdstA05(3) + zz*0, zz); 
hold on, plot(tt,-zz,'r'), legend('Data', 'Model')

subplot(2,4,7), plot([Wells.zo(4),Wells.z1000(4),Wells.z2000(4),Wells.z2500(4)],...
    [0 1000 2000 2500],'b'),set(gca,'Ydir','reverse'), title(Wells.Var3(5)),%ylim([-1000 1000]),xlim([10 150])
hold on, plot([Wells.zo(5),Wells.z1000(5),Wells.z2000(5),Wells.z2500(5)],[0 1000 2000 2500],'b'),
tt = FF(Wells.profdstA05(5) + zz*0, zz); 
hold on, plot(tt,-zz,'r'), 
plot(DemWelTVal.well4(:,1),DemWelTVal.well4(:,2)-DemWel.Elevation(4),'b'),
plot(DemWelTVal.well5(:,1),DemWelTVal.well5(:,2)-DemWel.Elevation(5),'b')