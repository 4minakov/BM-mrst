function obj = temperatureMismatch(model, states, varargin)

    loadData = @(name) load(fullfile(mrstPath('baia-mare'), name));
    
    data1 = loadData('Wells1.mat'); Wells = data1.Wells;
    data2 = loadData('Wells2.mat'); DemWel = data2.DemWel;

    z = -[0 1000 2000 2500]';
    x = repmat(DemWel.profdstA05(5), numel(z), 1);  

    x = [x, z];
    cells = findEnclosingCell(model.G, x);
    
    Tobs = [Wells.zo(5), Wells.z1000(5), Wells.z2000(5), Wells.z2500(5)]';
    
    Tsim = convertToCelcius(states{end}.T(cells));
    
    obj = {norm(Tobs - Tsim, 2)./273.15};
    
end