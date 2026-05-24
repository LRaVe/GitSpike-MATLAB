%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function [trials,labels] = build_trials(spikes,neurons)
    
    [N,S,R] = size(spikes);
    
    T = S*R;
    
    trials = cell(1,T);
    
    labels = zeros(1,T);
    
    idx = 0;
    
    for s = 1:S
    
        for r = 1:R
    
            idx = idx + 1;
    
            %% pooled train
    
            pooled = [];
    
            for n = neurons
    
                pooled = [pooled spikes{n,s,r}];
    
            end
    
            pooled = sort(pooled);
    
            trials{idx} = pooled;
    
            labels(idx) = s;
    
        end
    end
end


