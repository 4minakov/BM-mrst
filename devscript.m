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
% set up 'matrix' for parameter options for easier editing. The specific
% limits set for the various parameters influences the tuning/optimization 
% procedure to a large extent

setup0 = setup;
setup0.model.rock.permVal = setup.model.rock.perm;
setup0.model.rock.perm    = @(p,T) setup0.model.rock.permVal + p.*0;
setup0.model = setup0.model.setupOperators();

G = setup0.model.G;
pcell = G.cells.tag.crus_cell.*1 + G.cells.tag.sed1_cell.*2 + G.cells.tag.sed2_cell.*3;
pface = max(pcell(setup0.model.operators.N), [], 2);

pv = nan(max(pcell),1);
pv(pcell) = setup0.model.operators.pv;

parameters = [];
parameters = addParameter(parameters, setup0, ...
    'name', 'permeability', ...
    'location', {'rock', 'perm'}, ...
    'belongsTo', 'model', ...
    'setfun', @(varargin) setPermeability(varargin{:}), ...
    'scaling', 'log', ...
    'relativeLimits', [1e-4, 1e4], ...
    'lumping', pcell ...
);
parameters = addParameter(parameters, setup0, ...
    'name', 'porosity', ...
    'location', {'rock', 'poro'}, ...
    'belongsTo', 'model', ...
    'setfun', @(varargin) setPorosity(varargin{:}), ...
    'scaling', 'linear', ...
    'boxLims', [0.01, 0.8], ...
    'lumping', pcell ...
);
parameters = addParameter(parameters, setup0, ...
    'name', 'radiogenicHeatFluxDensity', ...
    'location', {'radiogenicHeatFluxDensity'}, ...
    'belongsTo', 'model', ...
    'scaling', 'linear', ...
    'boxLims', [0.001, 100]*micro*watt/meter^3, ...
    'lumping', pcell ...
);

pvec = getScaledParameterVector(setup0, parameters);

%% Model calibration Quasi-Newton
% make handle
loadData = @(name) load(fullfile(mrstPath('baia-mare'), name));
data     = loadData('Wells1.mat'    ); Wells1 = data.Wells;
data     = loadData('Wells2.mat'    ); DemWel = data.DemWel;
data     = loadData('DemWell1-6.mat'); DemWelTVal = data;

include = {'X3 -Cavnic', '125-VRO', '126-ROT'};
keep = ismember(Wells1.Var3, include);

Wells1 = Wells1(keep,:);

Tobs = [Wells1.zo, Wells1.z1000, Wells1.z2000, Wells1.z2500];
z0 = Wells1.Elevation;
z  = z0 -[0 1000 2000 2500];

x = repmat(Wells1.profdstA05, 1, 4);

Tobs = reshape(Tobs', [], 1);
x = reshape(x', [], 1);
z = reshape(z', [], 1);
x = [x,z];
cells = findEnclosingCell(setup.model.G, x);

inside = cells > 0;
cells = cells(inside);
Tobs = Tobs(inside);

%%
close all
mismatchFn = @(model, states, schedule, states_ref, compDer, tstep, state) ...
    temperatureMismatch(model, states, schedule, Tobs, cells, ...
                   'computePartials', compDer, 'tstep', tstep, ...
                   'state', state, 'from_states', false);
objh = @(p) evaluateMatch(p, mismatchFn, setup0, parameters, statesRef);
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 30 iterations
[objVal, pvecOpt, history] = unitBoxBFGS(pvec, objh, 'objChangeTol', 1e-5, ...
                                  'maxIt', 20, 'logPlot', true);

% %% Model calibration Levenberg-Marquard (using full Jacobian)
% mismatchFn2 = @(model, states, schedule, states_ref, compDer, tstep, state) ...
%     temperatureMismatch(model, states, schedule, Tobs, cells, ...
%         'computePartials', compDer, 'tstep', tstep,...
%         'state', state, 'from_states', false, 'mismatchSum', false);
% objh2 = @(p) evaluateMatchSummands(p, mismatchFn2, setup0, parameters, statesRef);
% % The calibration can be improved by taking a large number of iterations,
% % but here we set a limit of 30 iterations
% [v2, p_opt2, history2] = unitBoxLM(pvec, objh2, 'maxIt', 20);

%%
setupOpt = updateSetupFromScaledParameters(setup, parameters, pvecOpt); 
problemOpt = setupOpt.getPackedSimulationProblem();
simulatePackedProblem(problemOpt, 'restartStep', 1);
                              
%%
[~, states, reports] = getPackedSimulatorOutput(problemOpt);
setupOpt.plot(states);

%%
FF = scatteredInterpolant;
FF.Points = setup.model.G.cells.centroids; 
FF.ExtrapolationMethod='none';
FF.Values = convertToCelcius(states{end}.T);

%%
DemWelTVal = loadData('DemWell1-6.mat');
data1 = loadData('Wells1.mat'); Wells1 = data1.Wells;
data2 = loadData('Wells2.mat'); DemWel = data2.DemWel;

figure,  hold on,

plot(Wells1.profdstA05,Wells1.Elevation,'^','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10)
text(Wells1.profdstA05,Wells1.Elevation,Wells1.Var3)

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

subplot(2,4,3), plot([Wells1.zo(2),Wells1.z1000(2),Wells1.z2000(2),Wells1.z2500(2)],...
    [0 1000 2000 2500],'b'),set(gca,'Ydir','reverse'), title(Wells1.Var3(2)),%ylim([-1000 1000]),xlim([10 150])
tt = FF(Wells1.profdstA05(2) + zz*0, zz); 
hold on, plot(tt,-zz,'r'), legend('Data', 'Model')

subplot(2,4,4), plot([Wells1.zo(1),Wells1.z1000(1),Wells1.z2000(1),Wells1.z2500(1)],...
    [0 1000 2000 2500],'b'),set(gca,'Ydir','reverse'), title(Wells1.Var3(1)),%ylim([-1000 1000]),xlim([10 150])
hold on, plot(DemWelTVal.well6(:,1),DemWelTVal.well6(:,2)-DemWel.Elevation(6),'b'),
tt = FF(Wells1.profdstA05(1) + zz*0, zz); 
hold on, plot(tt,-zz,'r'), legend('Data', 'Model')

set(gca,'Ydir','reverse'), title(DemWel.Site(6)),%ylim([-1000 1000]),xlim([10 150])
subplot(2,4,5), plot([Wells1.zo(6),Wells1.z1000(6),Wells1.z2000(6),Wells1.z2500(6)],...
    [0 1000 2000 2500],'b'),set(gca,'Ydir','reverse'), title(Wells1.Var3(6)),%ylim([-1000 1000]),xlim([10 150])
tt = FF(Wells1.profdstA05(6) + zz*0, zz); 
hold on, plot(tt,-zz,'r'), legend('Data', 'Model')

subplot(2,4,6), plot([Wells1.zo(3),Wells1.z1000(3),Wells1.z2000(3),Wells1.z2500(3)],...
    [0 1000 2000 2500],'b'),set(gca,'Ydir','reverse'), title(Wells1.Var3(3)),%ylim([-1000 1000]),xlim([10 150])
tt = FF(Wells1.profdstA05(3) + zz*0, zz); 
hold on, plot(tt,-zz,'r'), legend('Data', 'Model')

subplot(2,4,7), plot([Wells1.zo(4),Wells1.z1000(4),Wells1.z2000(4),Wells1.z2500(4)],...
    [0 1000 2000 2500],'b'),set(gca,'Ydir','reverse'), title(Wells1.Var3(5)),%ylim([-1000 1000]),xlim([10 150])
hold on, plot([Wells1.zo(5),Wells1.z1000(5),Wells1.z2000(5),Wells1.z2500(5)],[0 1000 2000 2500],'b'),
tt = FF(Wells1.profdstA05(5) + zz*0, zz); 
hold on, plot(tt,-zz,'r'), 
plot(DemWelTVal.well4(:,1),DemWelTVal.well4(:,2)-DemWel.Elevation(4),'b'),
plot(DemWelTVal.well5(:,1),DemWelTVal.well5(:,2)-DemWel.Elevation(5),'b')

%%
function s = setPorosity(s, varargin)
    poro = varargin{end};
    s.rock.poro = poro;
    s.operators.pv = s.operators.vol.*poro;
end

function s = setPermeability(s, varargin)
    perm = varargin{end};
    s.rock.permVal = perm;
    s.rock.perm    = @(p,T) perm + p.*0;
end
