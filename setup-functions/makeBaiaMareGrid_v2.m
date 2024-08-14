clear all
close all
%% interpretation from geological transection

geol = readtable('geolCoordinates.txt');

img = imread('05_A-5.jpeg');

xleft  = 4804;
xright = 7250;
ytop   = 1488;
ybottom = 2109;

dist = distance(geol.Lat(2),geol.Lon(2),geol.Lat(3),geol.Lon(3),referenceEllipsoid('WGS84'));

figure('WindowState','maximized'), 
imagesc([0 dist]*1e-3, [-2,14], img(ytop:ybottom,xleft:xright,:)), hold on
xlabel('Distance [km]'), ylabel('Depth [km]')
set(gca,'Ydir','reverse','FontSize',14)

flist = dir('horizons/*.mat');
for i = 1:size({flist.name},2)
    load(['horizons/',flist(i).name]);
    fname = who(matfile(['horizons/',flist(i).name]));
    s = eval(fname{:});
    plot(s(:,1)*1e-3,s(:,2)*1e-3,'k','LineWidth',2)
end

%% grid setup

bottom = [0 7e3; topography(end,1) 7e3;];
ii = topography(:,1) > 0e3;
topography = topography(ii,:);
topography = [0e3, topography(1,2); topography];
outbnd = [topography; flipud(bottom); topography(1,:)];

% List of horizons
% ------------------
% lavaFlow1; % ok
% lavaFlow2; % ok 
% lavaFlow3; % ok
% intrusion1; % ok
% intrusion2; % failed
% intrusion3; % failed
% intrusion4; % ok
% baseMiocene1; % ok
% baseMiocene2; % failed
% basement; % ok
% sedimJK; % failed
% sedimTK; % failed
% sedimK2; % failed
% fault1; % ok
% fault2; % ok
% fault3; % ok
% fault4; % failed
% thrust1; % failed
% thrust2; % failed
% thrust3; % ok
% thrust4; % ok
% thrust5; % failed
% -------------------

surf = {baseMiocene1, basement, ...
    lavaFlow1, lavaFlow2, lavaFlow3, ...
    intrusion1, intrusion4, ...
    };

faults = {fault1, fault2, fault3, ...
    thrust3, thrust4};

rng(20240814);
n = 20; % Approximate number of cells in x-direction
L = max(outbnd)-min(outbnd);
dx = max(L)/n;
G = pebiGrid2D(dx, L, ...
    'polybdr'        , outbnd    , ...
    'FCFactor'       , 0.1       , ...
    'FCRefinement'   , true      , ...
    'FCEps'          , 0.1*max(L), ...
    'interpolateFC'  , true      , ...
    'CCFactor'       , 0.1       , ...
    'protLayer'      , true      , ...
    'protD'          , {@(x) dx*0.075}, ...
    'faceConstraints', surf      , ...
    'cellConstraints', faults     ...
    );
G = removeShortEdges(G, 1); % The grid may contain very short edges.

G = computeGeometry(G);     % Compute grid geometry
%%
figure
plotGrid(G,'faceColor','none'); axis equal tight, box on, hold on     
plot(outbnd(:,1),outbnd(:,2),'r','LineWidth',2),
xlabel('Distance [km]'), ylabel('Depth [km]')
set(gca,'Ydir','reverse','FontSize',14)
