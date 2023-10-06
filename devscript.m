mrstModule add baia-mare
mrstModule add test-suite
mrstModule add ad-core ad-props compositional
mrstModule add geothermal
mrstModule add upr

clc; close all

%% Set up test case
setup = TestCase('baia_mare_geothermal');
setup.plot();

%%
problem = setup.getPackedSimulationProblem();
simulatePackedProblem(problem);

%%
[~, states, reports] = getPackedSimulatorOutput(problem);
setup.plot(states);