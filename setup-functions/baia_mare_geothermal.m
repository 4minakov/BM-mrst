function setup = baia_mare_geothermal(varargin)
%Setup function for the Baia Mare basin

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    % Step 1: Test case description and options
    %---------------------------------------------------------------------%
    description = '';
    options = struct( ...
        'gridOpts', {{}}, ...
        'readGridFromDisk', true ... 
    );
    % Process optinal input arguments
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
                                        options, description, varargin{:});
    if ~fullSetup, return; end
    %---------------------------------------------------------------------%
    
    % Step 2: Define any module dependencies for the test case and set up
    %---------------------------------------------------------------------%
    % Module dependencies
    require ad-core ad-props geothermal compositional upr
    dataPath = fullfile(mrstPath('baia-mare'), 'data');
    % Set gravity in the neagtive y-direction
    gravity reset on
    gravity([0, -9.81]);
    
    % Get or construct grid
    %---------------------------------------------------------------------%
    G = [];
    if options.readGridFromDisk && exist(fullfile(dataPath, 'grid.mat'), 'file')
        data = load(fullfile(dataPath, 'grid.mat'));
        G = data.G;
    else
        G = makeBaiaMareGrid(options.gridOpts{:});
        save(fullfile(dataPath, 'grid.mat'), 'G');
    end
    %---------------------------------------------------------------------%
    
    % Make fluid
    %---------------------------------------------------------------------%
    rhoWS = 1000*kilogram/meter^3;
    fluid = initSimpleADIFluid( ...
        'mu'    , 1.0e-3, ...
        'rho'   , rhoWS , ...
        'phases', 'W'     ...
    );
    fluid = addThermalFluidProps(fluid, ... % Original fluid
        'Cp'     , 4.2e3*joule/(Kelvin*kilogram), ... % Specific heat capacity
        'lambdaF', 0.6*watt/(meter*Kelvin)      , ... % Thermal conductivity
        'cT'     , 2e-4*1/Kelvin                , ... % Thermal expansion coefficient
        'useEOS' , false                          ... % Use equation of state
    ); 
    %---------------------------------------------------------------------%
    
    % Make rock
    %---------------------------------------------------------------------%
    tag = G.cells.tag;
    % Set permeability
    perm                 = repmat(1e-17*meter^2, G.cells.num, 1);
    perm(tag.sediment1)  = 1e-16*meter^2; 
    perm(tag.sediment2)  = 1e-14*meter^2; 
    perm(tag.fault)      = 1e-12*meter^2; 
    perm(tag.intrusion1) = 1e-18*meter^2;
    perm(tag.intrusion2) = 1e-18*meter^2;
    % Set porosity
    poro                 = repmat(0.05, G.cells.num, 1);
    poro(tag.sediment1)  = 0.2;
    poro(tag.sediment2)  = 0.2;
    poro(tag.crust)      = 0.07; 
    poro(tag.fault)      = 0.3; 
    poro(tag.intrusion1) = 0.05;
    poro(tag.intrusion2) = 0.05;
    % Make rock
    rock = makeRock(G, perm, poro);
    
    % Specific heat capacity
    CpR = 1e3*joule/(Kelvin*kilogram); 
    % Thermal conductivity
    lambdaR = ones(G.cells.num,1)*1.5*watt/(meter*Kelvin);
    lambdaR(tag.sediment2)  = 2*watt/(meter*Kelvin); 
    lambdaR(tag.crust)      = 3*watt/(meter*Kelvin); 
    lambdaR(tag.intrusion1) = 3*watt/(meter*Kelvin); 
    lambdaR(tag.intrusion2) = 3*watt/(meter*Kelvin); 
    % Density
    rhoR                 = repmat(2600*kilogram/meter^3, G.cells.num, 1);
    rhoR(tag.sediment1)  = 2000*kilogram/meter^3; 
    rhoR(tag.sediment2)  = 2200*kilogram/meter^3; 
    rhoR(tag.intrusion1) = 2700*kilogram/meter^3;
    rhoR(tag.intrusion2) = 2700*kilogram/meter^3;
    % Set thermal properties to rock
    rock = addThermalRockProps( ...
        rock              , ... % Original rock
        'CpR'    , CpR    , ... % Specific heat capacity
        'lambdaR', lambdaR, ... % Thermal conductivity
        'rhoR'   , rhoR   , ... % Rock density
        'tau'    , 1        ... % Tortuosity
    );
    %---------------------------------------------------------------------%
    
    % Make model
    %---------------------------------------------------------------------%
    % Model limits
    K0   = 273.15*Kelvin; % Zero degrees Celcius
    pMax = 2000*barsa;    % Maximum pressure
    TMin = K0;            % Minimum temperature 
    TMax = K0 + 1300;     % Maximum temperature 
    % Make model
    model = GeothermalModel(G, rock, fluid);
    % Set limits
    model.minimumTemperature = TMin;
    model.maximumTemperature = TMax;
    model.maximumPressure    = pMax;
    % Set radiogenic heat flux density
    q = repmat(10*micro*watt/meter^3, G.cells.num, 1);
    q(tag.sediment1)  = 1*micro*watt/meter^3;
    q(tag.sediment2)  = 1*micro*watt/meter^3;
    q(tag.fault)      = 1*micro*watt/meter^3;
    q(tag.intrusion1) = 0.1*micro*watt/meter^3;
    q(tag.intrusion2) = 0.1*micro*watt/meter^3;
    % Set to model
    model.radiogenicHeatFluxDensity = q;
    %---------------------------------------------------------------------%
    
    % Set up initial conditions
    %---------------------------------------------------------------------%
    % Surface and intrusion temperatures
    TSurf = convertFromCelcius(20);
    TIntr = 600*Kelvin;
    % Set hydrostatic thermal equilibrium
    tag = G.faces.tag;
    north = find(tag.north);
    northX = [0; G.faces.centroids(north,1);
              max(G.nodes.coords(:,1))];
    northY = [G.faces.centroids(north(1),2); 
              G.faces.centroids(north,2)   ;
              G.faces.centroids(north(end),2)];
    cellX     = G.cells.centroids(:,1);
    cellYSurf = interp1(northX, northY,cellX);
    dz        = cellYSurf-G.cells.centroids(:,2);

    dpdz  = rhoWS*norm(model.gravity); % Hydrostatic gradient
    p0    = dpdz*dz;
    dTdz  = 45*Kelvin/(kilo*meter);    % Geothermal gradient
    T0Lin = TSurf + dTdz*dz;
    T0    = T0Lin;
    % [T1_ra,T1_ra_adv] = bm_geotherm(dz);
    % T0 = T1_ra_adv+K0;
    % T0(T0<10)=10;
    % T0 = T0+273;
    T0(G.cells.tag.intrusion1) = TIntr + K0;
    T0(G.cells.tag.intrusion2) = TIntr + K0;
    state0   = initResSol(G, p0, 1);        % Pressure and saturation 
    state0.T = T0;                          % Temperature
    %---------------------------------------------------------------------%
    
    % Set up boundary conditions
    %---------------------------------------------------------------------%
    tag = G.faces.tag;
    nn = nnz(tag.north);
    ns = nnz(tag.south);
    % Surface pressure
    pSurf = zeros(nn,1); % @@ Why zero for the last 2/3?
    pSurf(1:ceil(nn/3)) = 1*barsa;
    % Basal heat flux
    q = 60e-3*joule/second/meter^2;
    Q = q.*G.faces.areas(tag.south);
    % Set mass BCs
    bc = [];
    bc = addBC(bc, find(tag.north), 'pressure', pSurf); 
    bc = addBC(bc, find(tag.south), 'flux'    , 0    );
    bc = addBC(bc, find(tag.east ), 'flux'    , 0    );
    bc = addBC(bc, find(tag.west ), 'flux'    , 0    );
    % Set energy BCs
    [T, H] = deal(nan(numel(bc.face),1));
    T(1:nn) = TSurf; 
    H(nn + (1:ns)) = Q; 
    H(nn+ns+1:end) = 0; % @@ This imposes insulated ends - is that intentional?
    bc = addThermalBCProps(bc, 'T', T, 'Hflux', H);
    %---------------------------------------------------------------------%
    
    % Set up schedule
    %---------------------------------------------------------------------%
    timesteps = rampupTimesteps(3.5e6*year, 1e5*year);
    schedule  = simpleSchedule(timesteps, 'bc', bc);
    %---------------------------------------------------------------------%
    
    % Plotting
    %---------------------------------------------------------------------%
    plotOptions = {
        'Size'              , [1000, 350]             , ...
        'PlotBoxAspectRatio', [1.0000, 0.2657, 0.0001], ...
        'Box'               , true                    , ...
    };
    %---------------------------------------------------------------------%
    
    % Step 3: Pack test case setup
    %---------------------------------------------------------------------%
    setup = packTestCaseSetup(mfilename,              ...
                          'description', description, ...
                          'options'    , options    , ...
                          'state0'     , state0     , ...
                          'model'      , model      , ...
                          'schedule'   , schedule   , ...
                          'plotOptions', plotOptions);
    %---------------------------------------------------------------------%
    
end