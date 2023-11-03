function setup = baia_mare_geothermal_2(varargin)
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
    dataPath = fullfile(mrstPath('baia-mare'));
    % Set gravity in the neagtive y-direction
    gravity reset on
    gravity([0, -9.81]);
    
    % Get or construct grid
    %---------------------------------------------------------------------%
    G = [];
    if options.readGridFromDisk && exist(fullfile(dataPath, 'G_new.mat'), 'file')
        data = load(fullfile(dataPath, 'G_new.mat'));
        [G, G2D] = deal(data.G, data.G2D);
    else
%         G = makeBaiaMareGrid(options.gridOpts{:});
%         save(fullfile(dataPath, 'grid.mat'), 'G');
        error('Grid missing');
    end
    %---------------------------------------------------------------------%
    
    
    H = load(fullfile(dataPath, 'Profile5A.mat'));

    H1 = [];
    H1(:,1)  =  H.xnew;
    H1(:,2)  = -H.topography_new;
    H1(end,2)=  H1(end-1,2);
    %
    H2 = [];
    H2(:,1)  =  H.xnew;
    H2(:,2)  = -H.sediments1_new;
    H2(end,2)= H2(end-1,2);
    %
    H3 = [];
    H3(:,1)  =  H.xnew;
    H3(:,2)  = -H.sediments2_new;
    H3(end,2)=  H3(end-1,2);
    %
    Moho = [];
    Moho(:,1)  =  H.xnew;
    Moho(:,2)  = -H.basement_new;
    Moho(end,2)=  Moho(end-1,2);
    %
    Intr1 = H.intrusion1;
    Intr1 = [Intr1; Intr1(1,:)];
    Intr1(:,2)=-Intr1(:,2)-200;
    %
    Intr2 = H.intrusion2;
    Intr2 = [Intr2; Intr2(1,:)];
    Intr2(:,2)=-Intr2(:,2)-200;
    %
    Intr3 = H.intrusion3;
    Intr3 = [Intr3; Intr3(1,:)];
    Intr3(:,2)=-Intr3(:,2)-200;
    %
    Fault=[];
    Fault(:,2) = linspace(300, -6e3, 50)';
    Fault(:,1) = 30e3;

    %
    data = load(fullfile(dataPath, 'Intrusions_new.mat'));
    [Intr1, Intr2, Intr3] = deal(data.Intr1, data.Intr2, data.Intr3);
    outbnd = [H1(:,1), H1(:,2); flipud(Moho(:,1)), flipud(Moho(:,2)); H1(1,1), H1(1,2)];
    
    L(1)   = max(outbnd(:,1))-min(outbnd(:,1));
    L(2)   = max(outbnd(:,2))-min(outbnd(:,2));

    % Mapping for props and boundaries
    % units / faults and intrusions
    H2_x = flipud(H2(:,1));
    H2_y = flipud(H2(:,2));

    H3_x = flipud(H3(:,1));
    H3_y = flipud(H3(:,2));

    Moho_x =  flipud(Moho(:,1));
    Moho_y =  flipud(Moho(:,2));

    sed1_outbnd  = [H1(:,1), H1(:,2); H2_x, H2_y; H1(1,1), H1(1,2)];
    sed2_outbnd  = [H2(:,1), H2(:,2); H3_x, H3_y; H2(1,1), H2(1,2)];
    crus_outbnd  = [H3(:,1), H3(:,2); Moho_x, Moho_y; H3(1,1), H3(1,2)];
    fault       = find(G2D.cells.tag==1);

    id_f        = find(G.cells.centroids(fault)>10000); 
    fault       = fault(id_f);

    sed1_cell    = inpolygon(G.cells.centroids(:,1), ...
        G.cells.centroids(:,2), sed1_outbnd(:,1), sed1_outbnd(:,2));
    sed2_cell    = inpolygon(G.cells.centroids(:,1), ...
        G.cells.centroids(:,2), sed2_outbnd(:,1), sed2_outbnd(:,2));

    crus_cell   = inpolygon(G.cells.centroids(:,1), G.cells.centroids(:,2), crus_outbnd(:,1), crus_outbnd(:,2));
    intr1_cell  = inpolygon(G.cells.centroids(:,1), G.cells.centroids(:,2), Intr1(:,1), Intr1(:,2));
    intr2_cell  = inpolygon(G.cells.centroids(:,1), G.cells.centroids(:,2), Intr2(:,1), Intr2(:,2));
    intr3_cell  = inpolygon(G.cells.centroids(:,1), G.cells.centroids(:,2), Intr3(:,1), Intr3(:,2));

    G.cells.tag = struct();
    G.cells.tag.crus_cell  = crus_cell;
    G.cells.tag.sed1_cell  = sed1_cell;
    G.cells.tag.sed2_cell  = sed2_cell;
    G.cells.tag.intr1_cell = intr1_cell;
    G.cells.tag.intr2_cell = intr2_cell;
    G.cells.tag.intr3_cell = intr3_cell;

    % boundary 
    border = find(G.faces.neighbors(:,2)==0); 
    id_w   = find(G.faces.centroids(border,1)==0);
    West   = border(id_w); 
    id_e   = find(G.faces.centroids(border,1)==max(G.faces.centroids(:,1)));
    East   = border(id_e); 
    id_s   = find(G.faces.centroids(border,1)>0 & G.faces.centroids(border,1)<max(G.faces.centroids(:,1)) & ...
            G.faces.centroids(border,2)<-6000);
    South  = border(id_s);    
    id_n   = find(G.faces.centroids(border,1)>0 & G.faces.centroids(border,1)<max(G.faces.centroids(:,1)) & ...
            G.faces.centroids(border,2)>-6000);
    North  = border(id_n);    

    % check results 
    figure
    plotFaces(G2D, West, 'edgeColor', 'b');
    plotFaces(G2D, East, 'edgeColor', 'r');
    plotFaces(G2D, South, 'edgeColor', 'g');
    plotFaces(G2D, North, 'edgeColor', 'k');
    % 
    figure
    plotGrid(G2D,'faceColor','none'); axis equal tight, box on                
    plotGrid(G2D, sed1_cell , 'faceColor', 'y', 'faceAlpha', 0.5);  
    plotGrid(G2D, sed2_cell , 'faceColor', 'b', 'faceAlpha', 0.5);  
    plotGrid(G2D, crus_cell , 'faceColor', 'm');
    plotGrid(G2D, intr1_cell , 'faceColor', 'r', 'faceAlpha', 0.5);
    plotGrid(G2D, intr2_cell , 'faceColor', 'r', 'faceAlpha', 0.5);
    plotGrid(G2D, intr3_cell , 'faceColor', 'r', 'faceAlpha', 0.5);
    plotGrid(G2D, fault, 'faceColor', 'k','faceAlpha', 0.5);

    % set fluid structure properties
    rhoWS = 1000*kilogram/meter^3;
    fluid = initSimpleADIFluid('mu'    , 1.0e-3, ...
                               'rho'   , rhoWS , ...
                               'phases', 'W'   );
    fluid = addThermalFluidProps(fluid           , ... % Original fluid
                                 'Cp'     , 4.2e3*joule/(Kelvin*kilogram), ... % Specific heat capacity
                                 'lambdaF', 0.6*watt/(meter*Kelvin)  , ... % Thermal conductivity
                                 'useEOS' , false, ... % Use equation of state
                                 'cT'     , 2e-4*1/Kelvin );  % Thermal expansion coefficient

    % make rock structure
    perm = ones(G.cells.num, 1)*1e-17*meter^2;
    perm(sed1_cell) = 1e-16*meter^2; 
    perm(sed2_cell) = 1e-15*meter^2; 
    perm(fault) = 1e-13*meter^2; 
    perm(intr1_cell) = 1e-18*meter^2;
    perm(intr2_cell) = 1e-18*meter^2;
    perm(intr3_cell) = 1e-18*meter^2;
    poro = ones(G.cells.num,1)*0.05;
    poro(sed1_cell) = 0.2;
    poro(sed2_cell) = 0.2;
    poro(crus_cell) = 0.07; 
    poro(fault) = 0.3; 
    poro(intr1_cell) = 0.05;
    poro(intr2_cell) = 0.05;
    poro(intr3_cell) = 0.05;    
    rock = makeRock(G, perm, poro);

    CpR = ones(G.cells.num,1)*1e3*joule/(Kelvin*kilogram);
    Lambda = ones(G.cells.num,1)*1.5*watt/(meter*Kelvin);
    Lambda(sed2_cell)  = 2*watt/(meter*Kelvin); 
    Lambda(crus_cell)  = 3*watt/(meter*Kelvin); 
    Lambda(intr1_cell) = 3*watt/(meter*Kelvin); 
    Lambda(intr2_cell) = 3*watt/(meter*Kelvin); 
    Lambda(intr3_cell) = 3*watt/(meter*Kelvin); 
    RhoR  = ones(G.cells.num,1)*2600*kilogram/meter^3;
    RhoR(sed1_cell) = 2000*kilogram/meter^3; 
    RhoR(sed2_cell) = 2200*kilogram/meter^3; 
    RhoR(intr1_cell) = 2700*kilogram/meter^3;
    RhoR(intr2_cell) = 2700*kilogram/meter^3;
    RhoR(intr3_cell) = 2700*kilogram/meter^3;

    rock = addThermalRockProps(rock           , ... % Original rock
                               'CpR'    , CpR, ... % Specific heat capacity
                               'lambdaR', Lambda   , ... % Thermal conductivity
                               'rhoR'   , RhoR, ... % Rock density
                               'tau'    , 1   );    % Tortuosity
    %
    figure
    plotCellData(G, rock.lambdaR), axis equal, axis tight
    colorbar
    figure()
    plotCellData(G, rock.rhoR), axis equal, axis tight
    colorbar
    figure()
    plotCellData(G, log10(rock.perm)), axis equal, axis tight
    colorbar,alpha 0.5
    figure()
    plotCellData(G, rock.poro), axis equal, axis tight
    colorbar



    % Define the problem \
    % Define the gravity vector for a 2D x,y case using MRST convention

    % --- P and T frame
    K0   = 273.15*Kelvin; % Zero degrees Celcius
    pMin = 0*barsa;       % Minimum pressure
    pMax = 2000*barsa;    % Maximum pressure
    TMin = K0;            % Minimum temperature 
    TMax = K0 + 1300;      % Maximum temperature 

    gravity reset on
    gravity([0 -9.81]); % gravity is downward 
    model = GeothermalModel(G, rock, fluid);
    % Define limits of temperature validity for EOS
    model.minimumTemperature = TMin;
    model.maximumTemperature = TMax;
    model.maximumPressure    = pMax;
    %model.radiogenicHeatFluxDensity = 10*(micro*watt/meter^3);%radiogenic sources
    model.radiogenicHeatFluxDensity = ones(G.cells.num, 1)*2*(micro*watt/meter^3);%radiogenic sources
    model.radiogenicHeatFluxDensity(sed1_cell) = .1*(micro*watt/meter^3);%radiogenic sources 
    model.radiogenicHeatFluxDensity(sed2_cell) = .1*(micro*watt/meter^3);%radiogenic sources  
    model.radiogenicHeatFluxDensity(fault) = .1*(micro*watt/meter^3);%radiogenic sources 
    model.radiogenicHeatFluxDensity(intr1_cell) = .1*(micro*watt/meter^3);%radiogenic sources 
    model.radiogenicHeatFluxDensity(intr2_cell) = .1*(micro*watt/meter^3);%radiogenic sources 
    model.radiogenicHeatFluxDensity(intr3_cell) = .1*(micro*watt/meter^3);%radiogenic sources 
    %
    figure
    plotCellData(G, model.radiogenicHeatFluxDensity), axis equal, axis tight
    colorbar


    % Set up boundary conditions
    Tsurf  = K0+20;
    %psurf  = 1*barsa;
    psurf = zeros(numel(North),1);  
    psurf(1:ceil(numel(North)/3))=1*barsa;
    Tintr  = 600;
    heat_f = 60e-3*(joule/second/meter^2); 
    b_f    = heat_f.*G.faces.areas(South);
    %  BCs
    bc  = [];
    bc  = addBC(bc, North, 'pressure', psurf); 
    bc  = addBC(bc, South, 'flux', 0);
    bc  = addBC(bc, East, 'flux', 0);
    bc  = addBC(bc, West, 'flux', 0);
    T   = nan(numel(bc.face),1); 
    H   = nan(numel(bc.face),1);
    T(1:numel(North)) = Tsurf; 
    H(numel(North) +(1:numel(South))) = b_f; 
    H(numel(North) + numel(South) +1 : end) = 0;
    bc  = addThermalBCProps(bc, 'T', T, 'Hflux', H);


    % Setup initial conditions 
    North_x = [0; G.faces.centroids(North,1); max(L)];
    North_y = [G.faces.centroids(North(1),2); G.faces.centroids(North,2); G.faces.centroids(North(end),2)];
    cell_x  = G.cells.centroids(:,1);
    cell_y_surf = interp1(North_x, North_y,cell_x);
    dZ = cell_y_surf-G.cells.centroids(:,2);

    dpdz  = 1000*9.81; % rho*g  
    P0    = dpdz*dZ;
    dTdz  = 45*Kelvin/(kilo*meter); % Geothermal gradient
    T0_lin    = dTdz*dZ+ Tsurf;

    % @@ Need bm_geotherm!
%     [T1_ra,T1_ra_adv] = bm_geotherm(dZ);
%     T0 = T1_ra_adv+K0;
    
    T0 = T0_lin;
    
    %GeolPr = load('slice5A');
    %
    % FF =griddedInterpolant({GeolPr.prdepth1',GeolPr.prdist},GeolPr.slice5A1);
    % T0 = FF(-G.cells.centroids(:,2),G.cells.centroids(:,1)-110e3);
    T0(T0<10)=10;
    T0 = T0+273;
    T0(intr1_cell) = Tintr+K0;
    T0(intr2_cell) = Tintr+K0;
    T0(intr3_cell) = Tintr+K0;
    %
    state0   = initResSol(G, P0, 1);        % Pressure and saturation 
    state0.T = T0;                          % Temperature

    % Set up schedule
    timesteps = rampupTimesteps(3.5e6*year, 1e5*year); 
    schedule  = simpleSchedule(timesteps, 'bc', bc);
    
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