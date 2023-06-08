% Baia Mare 
clear all
close all

addpath('')
mrstModule add incomp linearsolvers mrst-gui ad-props  ad-core 
mrstModule add upr
mrstModule add ad-core ad-props geothermal compositional
mrstModule add mimetic incomp streamlines
mrstVerbose on 
%%
% load('HORIZ5A_2.mat')
% 
% HORIZ5A_2.H1(:,2) = -smooth(HORIZ5A_2.H1(:,2));
% HORIZ5A_2.H2(:,2) = -smooth(HORIZ5A_2.H2(:,2), 0.5);
% HORIZ5A_2.H3(:,2) = -smooth(HORIZ5A_2.H3(:,2));
% HORIZ5A_2.Moho(:,2) = -smooth(HORIZ5A_2.Moho(:,2)) + 20e3;
% 
% save HORIZ5A_3 HORIZ5A_2

load('HORIZ5A_3.mat')
%%

%
H1 = HORIZ5A_2.H1;
H1(:,1) = H1(:,1) - 110e3;
%H1(:,2) = -smooth(H1(:,2));
ii=find(H1(:,1)>0e3);
H1 = H1(ii,:);
H1 = [0e3, H1(1,2); H1];

%
H2 = HORIZ5A_2.H2 * 1000;
H2(:,1) = H2(:,1) -110e3;
%H2(:,2) = -smooth(H2(:,2),0.5);
ii=find(H2(:,1)>0);
H2 = H2(ii,:);
H2(:,2)=H2(:,2)-500;

load Sedim1
H2 = Sedim1;


H2 = [0e3, H2(1,2); H2];
H2(end,1) = H1(end,1);
%
H3 = HORIZ5A_2.H3;
H3(:,1) = H3(:,1) - 110e3;
%H3(:,2) = -smooth(H3(:,2));
ii=find(H3(:,1)>0e3);
H3 = H3(ii,:);
H3(:,2)=H3(:,2)-1000;
H3 = [0e3, H3(1,2); H3];
H3(end-1:end,:)=[];
H3(end,1) = H1(end,1);
%
Moho = HORIZ5A_2.Moho;
Moho(:,1) = Moho(:,1) - 110e3;
%Moho(:,2) = -smooth(Moho(:,2))+20e3;
ii=find(Moho(:,1)>0e3);
Moho = Moho(ii,:);
Moho = [0e3, Moho(1,2); Moho];
Moho(end,1) = H1(end,1);
%
load Intrus1
Intr1 = Intrus1;
% Intr1 = HORIZ5A_2.intrus1.*1000;
% Intr1(:,1) =  Intr1(:,1) - 110e3;
% Intr1(:,2) = -Intr1(:,2);
% ind = find(Intr1(:,2) < -12e3 | Intr1(:,2) > -4e3);
% Intr1(ind,:)=[];
Intr1 = [Intr1; Intr1(1,1), Intr1(1,2)];

load Intrus2
Intr2 = Intrus2;
Intr2 = [Intr2; Intr2(1,1), Intr2(1,2)];
%Intr2 = HORIZ5A_2.intrus2.*1000;
% Intr2(:,1) =  Intr2(:,1) - 110e3;
% Intr2(:,2) = -Intr2(:,2) ;
% ind = find(Intr2(:,2) < -12e3 | Intr2(:,2) > -4e3);
% Intr2(ind,:)=[];
% Intr2 = [Intr2; Intr2(1,1), Intr2(1,2)];
% Intr2(end-1,1) = Intr2(end-2,1);
%
Fault = HORIZ5A_2.fault_zone;
Fault(:,2) = -Fault(:,2);
Fault = [Fault(1,:); Fault(4,:)];
fault_y = linspace(500, -10e3, 50)'; 
fault_x = ones(length(fault_y), 1)*Fault(1,1)-110e3;
Fault = [fault_x, fault_y]; 


figure()
plot(H1(:,1), H1(:,2),'o-')
hold on 
plot(H2(:,1), H2(:,2),'o-')
plot(H3(:,1), H3(:,2),'o-')
plot(Moho(:,1), Moho(:,2),'o-')
plot(Fault(:,1), Fault(:,2),'o-')
plot(Intr1(:,1), Intr1(:,2),'o-r')
plot(Intr2(:,1), Intr2(:,2),'o-b')

%% Define upr grid 

Moho_x = flip(Moho(:,1));
Moho_y = flip(Moho(:,2));
outbnd = [H1(:,1), H1(:,2); Moho_x, Moho_y; H1(1,1), H1(1,2)];
% H2(1,:)=[]; H2(end,:)=[];
% H3(1,:)=[]; H3(end,:)=[];
surf = {H2, H3, Intr1, Intr2};
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
%     'CCFactor', 0.1 ); 
%  
% G2D = removeShortEdges(G2D, 1); % The grid may contain very short edges.
%% compute grid
%G = computeGeometry(G2D);
%save G G G2D
load G
%
figure 
plotGrid(G); axis equal tight, box on                 
plotFaces(G, G.faces.tag, 'edgeColor', 'b', 'lineWidth', 1); 
plotGrid(G, G.cells.tag , 'faceColor', 'r');                 


%% Mapping for props and boundaries
% units / faults and intrusions
H2_x = flip(H2(:,1));
H2_y = flip(H2(:,2));

H3_x = flip(H3(:,1));
H3_y = flip(H3(:,2));

sed1_outbnd  = [H1(:,1), H1(:,2); H2_x, H2_y; H1(1,1), H1(1,2)];
sed2_outbnd  = [H2(:,1), H2(:,2); H3_x, H3_y; H2(1,1), H2(1,2)];
crus_outbnd = [H3(:,1), H3(:,2); Moho_x, Moho_y; H3(1,1), H3(1,2)];
fault       = find(G2D.cells.tag==1);
% id_f        = find(G.cells.centroids(fault,1)>20000); 
% fault       = fault(id_f);

sed1_cell    = inpolygon(G.cells.centroids(:,1), ...
    G.cells.centroids(:,2), sed1_outbnd(:,1), sed1_outbnd(:,2));
sed2_cell    = inpolygon(G.cells.centroids(:,1), ...
    G.cells.centroids(:,2), sed2_outbnd(:,1), sed2_outbnd(:,2));

crus_cell   = inpolygon(G.cells.centroids(:,1), G.cells.centroids(:,2), crus_outbnd(:,1), crus_outbnd(:,2));
intr1_cell  = inpolygon(G.cells.centroids(:,1), G.cells.centroids(:,2), Intr1(:,1), Intr1(:,2));
intr2_cell  = inpolygon(G.cells.centroids(:,1), G.cells.centroids(:,2), Intr2(:,1), Intr2(:,2));


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
figure()
plotFaces(G2D, West, 'edgeColor', 'b');
plotFaces(G2D, East, 'edgeColor', 'r');
plotFaces(G2D, South, 'edgeColor', 'g');
plotFaces(G2D, North, 'edgeColor', 'k');
%% 
figure() 
plotGrid(G2D,'faceColor','none'); axis equal tight, box on                
plotGrid(G2D, sed1_cell , 'faceColor', 'y', 'faceAlpha', 0.5);  
plotGrid(G2D, sed2_cell , 'faceColor', 'b', 'faceAlpha', 0.5);  
plotGrid(G2D, crus_cell , 'faceColor', 'm');
plotGrid(G2D, intr1_cell , 'faceColor', 'r', 'faceAlpha', 0.5);
plotGrid(G2D, intr2_cell , 'faceColor', 'r', 'faceAlpha', 0.5);
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
perm(sed2_cell) = 1e-14*meter^2; 
perm(fault) = 1e-12*meter^2; 
perm(intr1_cell) = 1e-18*meter^2;
perm(intr2_cell) = 1e-18*meter^2;

poro = ones(G.cells.num,1)*0.05;
poto(sed1_cell) = 0.2;
poto(sed2_cell) = 0.2;
poro(crus_cell) = 0.07; 
poro(fault) = 0.3; 
poro(intr1_cell) = 0.05;
poro(intr2_cell) = 0.05;
    
rock = makeRock(G, perm, poro);

CpR = ones(G.cells.num,1)*1e3*joule/(Kelvin*kilogram);
Lambda = ones(G.cells.num,1)*1.5*watt/(meter*Kelvin);
Lambda(sed2_cell)  = 2*watt/(meter*Kelvin); 
Lambda(crus_cell)  = 3*watt/(meter*Kelvin); 
Lambda(intr1_cell) = 3*watt/(meter*Kelvin); 
Lambda(intr2_cell) = 3*watt/(meter*Kelvin); 
RhoR  = ones(G.cells.num,1)*2600*kilogram/meter^3;
RhoR(sed1_cell) = 2000*kilogram/meter^3; 
RhoR(sed2_cell) = 2200*kilogram/meter^3; 
RhoR(intr1_cell) = 2700*kilogram/meter^3;
RhoR(intr2_cell) = 2700*kilogram/meter^3;

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
model.radiogenicHeatFluxDensity = 10*(micro*watt/meter^3);%radiogenic sources

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
        quiver(G.cells.centroids(1:3:end,1),G.cells.centroids(1:3:end,2),...
        velocity(1:3:end,1),velocity(1:3:end,2),2,'k')
        title(['dt=',num2str(dt/(365*24*3600*1e6)),' my'])
        %caxis([270 350 ])
        colorbar, axis tight, shading faceted
        colormap(parula(32))
        pause(1)
        xlabel('Distance (m)'), ylabel('Depth (m)')
        drawnow
        set(gcf,'Units', 'normalized', 'Position', [ 0.1 0.1 0.7 0.6])
        %print('-dpng','-r300',[outfolder,'res',num2str(i)])
    end
    
    dt = dt + timesteps(i);
    
    
end

%%
x1 = linspace(0,50e3,300);
y1 = linspace(-12000,2000,300);
[x2d,y2d]=meshgrid(x1,y1);
figure(3),clf
dt = 0;

FF = scatteredInterpolant;
FF.Points = G.cells.centroids; 
FF.ExtrapolationMethod='none';
FF.Values = state0.T;
T2d0 = FF(x2d,y2d);
for i=1:length(timesteps)
    dt = dt + timesteps(i);
    if mod(i,1)==0
        FF.Values = states{i}.T;
        T2d = FF(x2d,y2d);
        T2d_old = T2d0;
        contourf(x2d,y2d,T2d-273,12), shading interp
        title(['dt=',num2str(dt/(365*24*3600*1e6)),' m.y.'])
        hold on
        velocity = faceFlux2cellVelocity(G, states{i}.flux);
        quiver(G.cells.centroids(1:3:end,1),G.cells.centroids(1:3:end,2),...
        velocity(1:3:end,1),velocity(1:3:end,2),2,'k')
        colorbar
        colormap(parula(32))
        %caxis([270 350])
        pause(0.5)
        xlabel('Distance (m)'), ylabel('Depth (m)')
        drawnow
        set(gcf,'Units', 'normalized', 'Position', [ 0.1 0.1 0.8 0.5])
        %print('-dpng','-r300',['out',num2str(i)])
        
    end
      
end
%%
% x1 = linspace(0,50e3,50);
% y1 = linspace(-12000,2000,50);
% [x2d,y2d]=meshgrid(x1,y1);
% figure(4),clf
% dt = 0;
% 
% FF = scatteredInterpolant;
% FF.Points = G.cells.centroids; 
% FF.ExtrapolationMethod='none';
% FF.Values = state0.pressure;
% P2d0 = FF(x2d,y2d);
% for i=1:length(timesteps)
%     dt = dt + timesteps(i);
%     if mod(i,1)==0
%         FF.Values = states{i}.pressure;
%         P2d = FF(x2d,y2d);
%         P2d_old = P2d0;
%         contourf(x2d,y2d,P2d,12), shading interp, hold on
%         
%         title(['dt=',num2str(dt/(365*24*3600*1e6)),' m.y.'])
%         colorbar
%         colormap(parula(32))
%         %caxis([270 350])
%         pause(.7)
%         xlabel('Distance (m)'), ylabel('Depth (m)')
%         drawnow
% %         set(gcf,'Units', 'normalized', 'Position', [ 0.1 0.1 0.7 0.3])
%         %print('-dpng','-r300',['Pout',num2str(i)])
%         
%     end
%       
% end
%%

