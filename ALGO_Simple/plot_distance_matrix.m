%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function plot_distance_matrix(D,labels,titleStr)
    
    P = compute_discrimination_performance(D,labels);

    imagesc(D);
    
    colormap jet;
    
    colorbar;
    
    title(sprintf('%s   |   P = %.3f',titleStr,P), 'FontWeight', 'bold');    

    xlabel('Trials');
    ylabel('Trials');
    
    hold on;
    
    %% =====================================
    %% lignes séparation stimuli
    %% =====================================
    
    S = max(labels);
    
    R = sum(labels==1);
    
    for s = 1:S-1
    
        pos = s*R + 0.5;
    
        xline(pos,'w','LineWidth',2);
        yline(pos,'w','LineWidth',2);
    
    end

end


