%% Model tuning of Baia Mare

% Add modules
mrstModule add baia-mare
mrstModule add test-suite
mrstModule add ad-core ad-props compositional
mrstModule add geothermal
mrstModule add upr
mrstModule add optimization
mrstModule add mrst-gui
% Command window output
mrstVerbose on
% Load datasets
loadData = @(name) load(fullfile(mrstPath('baia-mare'), name));

%% Set up test case
setup = TestCase('baia_mare_geothermal_2', 'readGridFromDisk', true);
setup.plot();

%% Simulate base case
problem = setup.getPackedSimulationProblem();
simulatePackedProblem(problem, 'restartStep', nan);

%% Get base case results
[~, states] = getPackedSimulatorOutput(problem);
T0 = states{end}.T;
setup.plot(states);
%% Set up optimization problem
% Initial guess
setup0 = setup;
% Set permeability as function
setup0.model.rock.permVal = setup.model.rock.perm;
setup0.model.rock.perm    = @(p,T) setup0.model.rock.permVal + p.*0;
setup0.model              = setup0.model.setupOperators();
% Make partition vector defining regions with equal parameters. Adjust
% this to a suitalble partition
G = setup0.model.G;
pcell = G.cells.tag.crus_cell.*1 ...
      + G.cells.tag.sed1_cell.*2 ...
      + G.cells.tag.sed2_cell.*3;
setup.plot(pcell);
% Set up parameters for tuning
% boxLimits: Minimum/maximum allowed parameter value
parameters = [];

% Permeability
%-------------------------------------------------------------------------%
permLim = [0.01, 2000; 0.01, 2000; 0.01, 2000]*milli*darcy;
parameters = addParameter(parameters, setup0, ...
    'name'     , 'permeability'                          , ...
    'location' , {'rock', 'perm'}                        , ...
    'belongsTo', 'model'                                 , ...
    'setfun'   , @(varargin) setPermeability(varargin{:}), ...
    'scaling'  , 'log'                                   , ...
    'boxLims'  , permLim                                 , ...
    'lumping'  , pcell                                     ...
);
%-------------------------------------------------------------------------%

% Porosity
%-------------------------------------------------------------------------%
poroLim = [0.001, 0.8; 0.001, 0.8; 0.001, 0.8];
parameters = addParameter(parameters, setup0, ...
    'name'     , 'porosity'                          , ...
    'location' , {'rock', 'poro'}                    , ...
    'belongsTo', 'model'                             , ...
    'setfun'   , @(varargin) setPorosity(varargin{:}), ...
    'scaling'  , 'linear'                            , ...
    'boxLims'  , poroLim                             , ...
    'lumping'  , pcell                                 ...
);
%-------------------------------------------------------------------------%

% Radiogenic heat flux density
%-------------------------------------------------------------------------%
qLim = [0.001, 100; 0.001, 100; 0.001, 100]*micro*watt/meter^3;
parameters = addParameter(parameters, setup0, ...
    'name'     , 'radiogenicHeatFluxDensity'  , ...
    'location' , {'radiogenicHeatFluxDensity'}, ...
    'belongsTo', 'model'                      , ...
    'scaling'  , 'linear'                     , ...
    'boxLims'  , qLim                         , ...
    'lumping'  , pcell                          ...
);
%-------------------------------------------------------------------------%

% Get parameter vector
pvec = getScaledParameterVector(setup0, parameters);

%% Load reference data
% The observed temperature and the corresponding grid cells must be
% provided to the temperatureMismatch (see below)
data     = loadData('Wells1.mat'    ); Wells1 = data.Wells;
data     = loadData('Wells2.mat'    ); DemWel = data.DemWel;
data     = loadData('DemWell1-6.mat'); DemWelTVal = data;

include = {'X3 -Cavnic', '125-VRO', '126-ROT'};
keep = ismember(Wells1.Var3, include);

Wells1 = Wells1(keep,:);

Tobs = [Wells1.zo, Wells1.z1000, Wells1.z2000, Wells1.z2500];
Tobs0 = Tobs;
ztop = Wells1.Elevation;
z    = ztop -[0 1000 2000 2500];
z0   = z;

x = repmat(Wells1.profdstA05, 1, 4);

Tobs = reshape(Tobs', [], 1);
x    = reshape(x', [], 1);
z    = reshape(z', [], 1);
x    = [x,z];
cells = findEnclosingCell(setup.model.G, x);
% Only match againts values measured inside simulation grid
inside = cells > 0;

cells = cells(inside);
Tobs  = Tobs(inside);

%% Perform model Quasi-Newton model calibration
close all
% Set up mismatch function
mismatchFn = @(model, states, schedule, states_ref, compDer, tstep, state) ...
    temperatureMismatch(model, states, schedule, Tobs, cells, ...
    'computePartials', compDer, 'tstep', tstep, ...
        'state', state, 'from_states', false);
objh = @(p) evaluateMatch(p, mismatchFn, setup0, parameters, []);
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 20
[objVal, pvecOpt, history] = unitBoxBFGS(pvec, objh, 'objChangeTol', 1e-5, ...
                                  'maxIt', 20, 'logPlot', true);

%% Make new setup with tuned parameters and simulate
setupOpt = updateSetupFromScaledParameters(setup, parameters, pvecOpt); 
problemOpt = setupOpt.getPackedSimulationProblem();
simulatePackedProblem(problemOpt, 'restartStep', 1);
                              
%% Load results from tuned model
[~, statesOpt, reports] = getPackedSimulatorOutput(problemOpt);
TOpt = convertToCelcius(statesOpt{end}.T);

%% Plot match
close all
fg = setup.figure();
plotGrid(setup.model.G); setup.setAxisProperties(gca), hold on,

colors = lines(numel(include));
for i = 1:numel(include)
    
    color = colors(i,:);
    figure(fg),
    plot(Wells1.profdstA05(i), Wells1.Elevation(i), ...
        'v', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k', 'MarkerSize', 10)
    text(Wells1.profdstA05(i), Wells1.Elevation(i), Wells1.Var3{i})
    
    ft = figure(); hold on
    To = Tobs0(i,:);
    z  = z0(i,:);
    plot(To, z, 'o', 'Color', color, 'LineWidth', 2);

    x = Wells1.profdstA05(i);
    n = 25;
    t = linspace(0,1,n)';    
    x = interp1([0,1]', [x,z(1); x, z(end)], t);
    
    cells = findEnclosingCell(setup.model.G, x);
    cells = cells(cells > 0);
    cells = rlencode(cells,1);
    
    zg = setup.model.G.cells.centroids(cells,2);
    plot(TOpt(cells), zg, 'LineWidth', 2, 'Color', color);
    axis tight, grid on, box on
    legend({'Measured', 'Simulated'});    
    figure(fg), plotGrid(setup.model.G, cells, 'faceColor', color)
    
end

%% Setters for permeability and porosity
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
