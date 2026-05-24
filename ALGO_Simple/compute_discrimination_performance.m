%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function P = compute_discrimination_performance(D,labels)
    
    T = length(labels);
    
    intra = [];
    inter = [];
    
    for i = 1:T
    
        for j = i+1:T
    
            %% =====================================
            %% same stimulus
            %% =====================================
    
            if labels(i) == labels(j)
    
                intra(end+1) = D(i,j);
    
            %% =====================================
            %% different stimuli
            %% =====================================
    
            else
    
                inter(end+1) = D(i,j);
    
            end
        end
    end
    
    %% =====================================
    %% discrimination performance
    %% =====================================
    
    P = mean(inter) - mean(intra);
    
end


