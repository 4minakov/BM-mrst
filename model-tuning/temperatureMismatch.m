function obj = temperatureMismatch(model, states, schedule, Tobs, cells, varargin)
% Compute mismatch-function 

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
opt     = struct('ComputePartials',     false, ...
                 'tStep' ,              [], ...
                 'state',               [],...
                 'from_states',         true,...% can be false for generic models
                 'matchOnlyProducers',  false, ...
                 'mismatchSum',         true, ...
                 'accumulateWells',       [], ...
                 'accumulateTypes',       []);
             
opt     = merge_options(opt, varargin{:});

dts   = schedule.step.val;

tSteps = opt.tStep;
if isempty(tSteps) %do all
    numSteps = numel(dts);
    tSteps = (1:numSteps)';
else
    numSteps = 1;
end

obj = repmat({[]}, numSteps, 1);

for step = 1:numSteps
   
    stepNo = tSteps(step);
    if opt.ComputePartials
        if(opt.from_states)
            init=true;
            state = model.getStateAD(states{stepNo}, init);
        else
            state = opt.state;
        end
    else
        state = states{stepNo};
    end
    T = model.getProp(state, 'temperature');
    T = convertToCelcius(T(cells));
    
    o = ((T-Tobs)./Tobs).^2;
    
    if opt.mismatchSum, o = sum(o); end
        
    if (stepNo < numel(schedule.step.val)), o = o.*0; end
    
    obj{step} = o;
end

end
%--------------------------------------------------------------------------