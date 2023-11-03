function G = makeBaiaMareGrid(varargin)
% Make grid of the Baia Mare basin
% @@ TODO: use padding for the fault cells

    % Get horizon, fault and intrusion data
    %---------------------------------------------------------------------%
    loadData = @(name) load(fullfile(mrstPath('baia-mare'), name)); %#ok

    horizons = loadData('HORIZ5A_3.mat');
    horizons = horizons.HORIZ5A_2;
    
    % Horizon 1
    h1 = horizons.H1;
    h1(:,1) = h1(:,1) - 110e3;
    ii = h1(:,1) > 0e3;
    h1 = h1(ii,:);
    h1 = [0e3, h1(1,2); h1];
    
    % Horizon 2
    sedim1 = loadData('Sedim1');
    h2 = sedim1.Sedim1;
    h2 = [0e3, h2(1,2); h2];
    h2(end,1) = h1(end,1);
    
    % Horizon 3
    h3 = horizons.H3;
    h3(:,1) = h3(:,1) - 110e3;
    ii = h3(:,1) > 0e3;
    h3 = h3(ii,:);
    h3(:,2) = h3(:,2)-1000;
    h3 = [0e3, h3(1,2); h3];
    h3(end-1:end,:)=[];
    h3(end,1) = h1(end,1);
    
    % Horizon Moho(?)
    moho = horizons.Moho;
    moho(:,1) = moho(:,1) - 110e3;
    ii = moho(:,1)>0e3;
    moho = moho(ii,:);
    moho = [0e3, moho(1,2); moho];
    moho(end,1) = h1(end,1);
    
    % Intrusion 1
    data = loadData('Intrus1');
    i1 = data.Intrus1;
    i1 = [i1; i1(1,1), i1(1,2)];

    % Intrusion 2
    data = loadData('Intrus2');
    i2 = data.Intrus2;
    i2 = [i2; i2(1,1), i2(1,2)];
    
    % Fault
    fault = horizons.fault_zone;
    fault(:,2) = -fault(:,2);
    fault   = [fault(1,:); fault(4,:)];
    fault_y = linspace(500, -10e3, 50)'; 
    fault_x = ones(length(fault_y), 1)*fault(1,1)-110e3;
    fault   = [fault_x, fault_y]; 
    %---------------------------------------------------------------------%
    
    % Construct PEBI grid 
    %---------------------------------------------------------------------%
    outbnd = [h1; flipud(moho); h1(1,1), h1(1,2)];
    surf = {h2, h3, i1, i2};
    fault  = {fault};

    rng(20231006);
    n = 20; % Approximate number of cells in x-direction
    L = max(outbnd)-min(outbnd);
    dx = max(L)/n;
    G = pebiGrid2D(dx, L, ...
        'polybdr'        , outbnd    , ... 
        'faceConstraints', surf      , ...
        'FCFactor'       , 0.1       , ...
        'FCRefinement'   , true      , ...
        'FCEps'          , 0.1*max(L), ...
        'interpolateFC'  , true      , ...
        'cellConstraints', fault     , ...
        'CCFactor'       , 0.1       , ...
        'protLayer'      , true      , ...
        'protD'          , {@(x) dx*0.075} ...
    ); 
    G = removeShortEdges(G, 1); % The grid may contain very short edges.
    G = computeGeometry(G);     % Compute grid geometry
    %---------------------------------------------------------------------%
    
    % Set cell tags
    %---------------------------------------------------------------------%
    % Find cells for each geoelogical region
    s1Bdr      = [h1; flipud(h2); h1(1,:)];
    s2Bdr      = [h2; flipud(h3); h2(1,:)];
    crustBdr   = [h3; flipud(moho); h3(1,:)];
    xc         = G.cells.centroids;
    s1Cells    = inpolygon(xc(:,1), xc(:,2), s1Bdr(:,1), s1Bdr(:,2));
    s2Cells    = inpolygon(xc(:,1), xc(:,2), s2Bdr(:,1), s2Bdr(:,2));
    crustCells = inpolygon(xc(:,1), xc(:,2), crustBdr(:,1), crustBdr(:,2));
    i1Cells    = inpolygon(xc(:,1), xc(:,2), i1(:,1), i1(:,2));
    i2Cells    = inpolygon(xc(:,1), xc(:,2), i2(:,1), i2(:,2));
    faultCells = G.cells.tag==1 & G.cells.centroids(:,1) > 10000; 
    % Make tag struct
    tag = struct();
    % Sediment 1
    tag.sediment1          = false(G.cells.num, 1);
    tag.sediment1(s1Cells) = true;
    % Sediment 2
    tag.sediment2          = false(G.cells.num, 1);
    tag.sediment2(s2Cells) = true;
    % Crust
    tag.crust             = false(G.cells.num, 1);
    tag.crust(crustCells) = true;
    % Intrusion 1
    tag.intrusion1        = false(G.cells.num, 1);
    tag.intrusion1(i1Cells) = true;
    % Intrusion 2
    tag.intrusion2        = false(G.cells.num, 1);
    tag.intrusion2(i2Cells) = true;
    % Fault
    tag.fault             = false(G.cells.num, 1);
    tag.fault(faultCells) = true;
    % Set tags to grid
    G.cells.tag = tag;
    %---------------------------------------------------------------------%
    
    % Set face tags
    %---------------------------------------------------------------------%
    tag = struct();
    border    = any(G.faces.neighbors == 0, 2); 
    tag.west  = border & G.faces.centroids(:,1) == 0;
    tag.east  = border & G.faces.centroids(:,1) == max(G.faces.centroids(:,1));
    tag.south = border & G.faces.centroids(:,1) > 0 ...
        & G.faces.centroids(:,1) < max(G.faces.centroids(:,1)) ...
        & G.faces.centroids(:,2) < -6000;
    tag.north = border & G.faces.centroids(:,1) > 0 ...
        & G.faces.centroids(:,1) < max(G.faces.centroids(:,1)) ...
        & G.faces.centroids(:,2) > -6000;
    G.faces.tag = tag;
    %---------------------------------------------------------------------%
    
end