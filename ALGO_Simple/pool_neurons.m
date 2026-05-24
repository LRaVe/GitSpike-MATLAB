%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function pooled = pool_neurons(spikes,neurons,s,r)
    
    pooled = [];
    
    for n = neurons
    
        pooled = [pooled spikes{n,s,r}];
    
    end
    
    pooled = sort(pooled);

end

