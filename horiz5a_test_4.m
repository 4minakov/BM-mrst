% Baia Mare 
clear all
close all


mrstPath reregister geothermal 'C:\Users\alexamin\Dropbox (UiO)\BaiaMare\fem3d\mrst-2022a\modules\mrst-geothermal'


addpath('')
mrstModule add incomp linearsolvers mrst-gui ad-props  ad-core 
mrstModule add upr
mrstModule add ad-core ad-props geothermal compositional
mrstModule add mimetic incomp streamlines
mrstVerbose on 

%% 
H = load('C:\Users\alexamin\Dropbox (UiO)\BaiaMare\fem3d\mrst-2022a\data\Profile5A')

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
%%
figure
hold on 
plot(H1(:,1), H1(:,2),'o-')
plot(H2(:,1), H2(:,2),'o-')
plot(H3(:,1), H3(:,2),'o-')
plot(Moho(:,1), Moho(:,2),'o-')
plot(Fault(:,1), Fault(:,2),'o-')
plot(Intr1(:,1), Intr1(:,2),'o-r')
plot(Intr2(:,1), Intr2(:,2),'o-b')
plot(Intr3(:,1), Intr3(:,2),'o-b')



%% Define upr grid 
load Intrusions_new
outbnd = [H1(:,1), H1(:,2); flipud(Moho(:,1)), flipud(Moho(:,2)); H1(1,1), H1(1,2)];
% H2(1,:)=[]; H2(end,:)=[];
% H3(1,:)=[]; H3(end,:)=[];
surf = {H2, H3, Intr1, Intr2, Intr3};
fault  = {Fault};

%% Generate Pebi grid
rng('default')
n   = 20; % Approximate number of cells in x-direction
L(1)   = max(outbnd(:,1))-min(outbnd(:,1));
L(2)   = max(outbnd(:,2))-min(outbnd(:,2));

% G2D = pebiGrid2D(max(L)/n, L           , ...
%      'polybdr'        , outbnd      , ... 
%      'faceConstraints', surf, ...
%     'FCFactor'       , 0.1 , ...
%     'FCRefinement'   , true ,      ...
%     'FCEps', 0.1*max(L), ...
%     'interpolateFC'  , true, ...
%     'cellConstraints', fault, ...
%       'CCFactor', 0.1 ); 
%  
% G2D = removeShortEdges(G2D, 1); % The grid may contain very short edges.
%% compute grid
%G = computeGeometry(G2D);
%save G_new G G2D
load G_new
%
figure 
plotGrid(G); axis equal tight, box on                 
plotFaces(G, G.faces.tag, 'edgeColor', 'b', 'lineWidth', 1); 
plotGrid(G, G.cells.tag , 'faceColor', 'r');                 


%% Mapping for props and boundaries
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


%% boundary 
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

%% check results 
figure
plotFaces(G2D, West, 'edgeColor', 'b');
plotFaces(G2D, East, 'edgeColor', 'r');
plotFaces(G2D, South, 'edgeColor', 'g');
plotFaces(G2D, North, 'edgeColor', 'k');
%% 
figure
plotGrid(G2D,'faceColor','none'); axis equal tight, box on                
plotGrid(G2D, sed1_cell , 'faceColor', 'y', 'faceAlpha', 0.5);  
plotGrid(G2D, sed2_cell , 'faceColor', 'b', 'faceAlpha', 0.5);  
plotGrid(G2D, crus_cell , 'faceColor', 'm');
plotGrid(G2D, intr1_cell , 'faceColor', 'r', 'faceAlpha', 0.5);
plotGrid(G2D, intr2_cell , 'faceColor', 'r', 'faceAlpha', 0.5);
plotGrid(G2D, intr3_cell , 'faceColor', 'r', 'faceAlpha', 0.5);
plotGrid(G2D, fault, 'faceColor', 'k','faceAlpha', 0.5);

%% set fluid structure properties
rhoWS = 1000*kilogram/meter^3;
fluid = initSimpleADIFluid('mu'    , 1.0e-3, ...
                           'rho'   , rhoWS , ...
                           'phases', 'W'   );
fluid = addThermalFluidProps(fluid           , ... % Original fluid
                             'Cp'     , 4.2e3*joule/(Kelvin*kilogram), ... % Specific heat capacity
                             'lambdaF', 0.6*watt/(meter*Kelvin)  , ... % Thermal conductivity
                             'useEOS' , false, ... % Use equation of state
                             'cT'     , 2e-4*1/Kelvin );  % Thermal expansion coefficient

%% make rock structure
perm = ones(G.cells.num, 1)*1e-17*meter^2;
perm(sed1_cell) = 1e-16*meter^2; 
perm(sed2_cell) = 1e-15*meter^2; 
perm(fault) = 1e-13*meter^2; 
perm(intr1_cell) = 1e-18*meter^2;
perm(intr2_cell) = 1e-18*meter^2;
perm(intr3_cell) = 1e-18*meter^2;
poro = ones(G.cells.num,1)*0.05;
poto(sed1_cell) = 0.2;
poto(sed2_cell) = 0.2;
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
%%
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



%% Define the problem \
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


%% Set up boundary conditions
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


%% Setup initial conditions 
North_x = [0; G.faces.centroids(North,1); max(L)];
North_y = [G.faces.centroids(North(1),2); G.faces.centroids(North,2); G.faces.centroids(North(end),2)];
cell_x  = G.cells.centroids(:,1);
cell_y_surf = interp1(North_x, North_y,cell_x);
dZ = cell_y_surf-G.cells.centroids(:,2);

dpdz  = 1000*9.81; % rho*g  
P0    = dpdz*dZ;
dTdz  = 45*Kelvin/(kilo*meter); % Geothermal gradient
T0_lin    = dTdz*dZ+ Tsurf;
[T1_ra,T1_ra_adv] = bm_geotherm(dZ);
T0 = T1_ra_adv+K0;
%GeolPr = load('slice5A');
%%
% FF =griddedInterpolant({GeolPr.prdepth1',GeolPr.prdist},GeolPr.slice5A1);
% T0 = FF(-G.cells.centroids(:,2),G.cells.centroids(:,1)-110e3);
T0(T0<10)=10;
T0 = T0+273;
T0(intr1_cell) = Tintr+K0;
T0(intr2_cell) = Tintr+K0;
T0(intr3_cell) = Tintr+K0;
%%
state0   = initResSol(G, P0, 1);        % Pressure and saturation 
state0.T = T0;                          % Temperature

%% Set up schedule
timesteps = rampupTimesteps(3.5e6*year, 1e5*year); 
schedule  = simpleSchedule(timesteps, 'bc', bc);
%% 
figure()
plotToolbar(G, state0, 'edgealpha', 0.2);
colorbar

%% Run simulation
[~, states,schedulereport] = simulateScheduleAD(state0, model, schedule);
%%

outfolder = 'C:\Users\alexamin\Dropbox (UiO)\BaiaMare\fem3d\mrst-2022a\marine\out\';

figure(2),clf
dt = 0;
for i=1:length(timesteps)
    if mod(i,2)==0
        states1 = states{i};
        states1.T = states1.T-273;
        plotToolbar(G, states1, 'edgealpha', 0.2, 'field','T');
         hold on
        velocity = faceFlux2cellVelocity(G, states{i}.flux);
        HFd = faceFlux2cellVelocity(G,states1.heatFluxCond);
%         HFa = faceFlux2cellVelocity(G,states1.heatFluxAdv);
        
        quiver(G.cells.centroids(1:3:end,1),G.cells.centroids(1:3:end,2),...
        velocity(1:3:end,1),velocity(1:3:end,2),4,'k')
        title(['dt=',num2str(dt/(365*24*3600*1e6)),' my'])
        %caxis([270 350 ])
        colorbar, axis tight, shading faceted
        colormap(parula(32))
        pause(.1)
        xlabel('Distance (m)'), ylabel('Depth (m)'),%ylim([-12000 2000])
        
        set(gcf,'Units', 'normalized', 'Position', [ 0.1 0.1 0.7 0.6])
        drawnow
        %print('-dpng','-r300',[outfolder,'res_noRA_',num2str(i)])
    end
    
    dt = dt + timesteps(i);
    
    
end

%%
figure,plotCellData(G, HFd(:,2)), caxis([30 200]*1e-3), colorbar
FF = scatteredInterpolant;
FF.Points = G.cells.centroids; 
FF.ExtrapolationMethod='none';
FF.Values = states1.T;
%%
DemWelTVal = load('C:\Users\alexamin\Dropbox (UiO)\BaiaMare\data\HeatFlow\DemWell1-6')
load 'C:\Users\alexamin\Dropbox (UiO)\BaiaMare\fem3d\mrst-2022a\data\Wells1'
load 'C:\Users\alexamin\Dropbox (UiO)\BaiaMare\fem3d\mrst-2022a\data\Wells2'
figure(2),  hold on,
plot(Wells.profdstA05,Wells.Elevation,'^','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10)
text(Wells.profdstA05,Wells.Elevation,Wells.Var3)

plot(DemWel.profdstA05,DemWel.Elevation,'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10)
text(DemWel.profdstA05,DemWel.Elevation,DemWel.Site), ylim([-9e3 2e3])
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
%%
