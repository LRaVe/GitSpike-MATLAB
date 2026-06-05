%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function P = compute_discrimination_performance(D,labels)
    
    T = length(labels);
    
    intra = zeros(1,2*T);
    inter = zeros(1,(T*(T-1)/2)-(2*T));
    idx_intra = 1;
    idx_inter = 1;

    for i = 1:T
    
        for j = i+1:T
    
            %% =====================================
            %% same stimulus
            %% =====================================
    
            if labels(i) == labels(j)

                intra(idx_intra) = D(i,j);
                idx_intra = idx_intra + 1;
    
            %% =====================================
            %% different stimuli
            %% =====================================
    
            else
    
                inter(idx_inter) = D(i,j);
                idx_inter = idx_inter + 1;
    
            end
        end
    end

    %% =====================================
    %% discrimination performance
    %% =====================================
    
    P = mean(inter) - mean(intra);
    
end


