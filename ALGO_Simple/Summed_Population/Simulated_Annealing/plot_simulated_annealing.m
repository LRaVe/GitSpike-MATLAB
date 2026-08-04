%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function plot_simulated_annealing(result,codingNeurons)
    
    figure('Color','w');
    
    tiledlayout(3,2,'TileSpacing','compact');
    
    %% ==========================================
    %% POPULATION SIZE
    %% ==========================================
    
    nexttile
    
    plot(result.history.size,'b','LineWidth',1.5);
    
    hold on
    
    yline(length(codingNeurons),'--r','Ground truth');
    
    xlabel('Iteration');
    ylabel('Population size');
    
    title('Population size');
    
    grid on
    
    
    %% ==========================================
    %% PERFORMANCE
    %% ==========================================
    
    nexttile
    
    plot(result.history.P,'k','LineWidth',1.5);
    
    hold on
    
    plot(result.history.bestP,'r','LineWidth',2);
    
    xlabel('Iteration');
    ylabel('P');
    
    legend({'Current','Best'},'Location','best');
    
    title('Discrimination performance');
    
    grid on
    
    
    %% ==========================================
    %% CODING / NON-CODING NEURONS
    %% ==========================================
    
    nexttile
    
    plot(result.history.nCoding,'g','LineWidth',1.5);
    
    hold on
    
    plot(result.history.nNonCoding,'r','LineWidth',1.5);
    
    yline(length(codingNeurons),'--g','All coding neurons');
        
    xlabel('Iteration');
    
    ylabel('Neuron count');
    
    title('Population composition');
    
    legend({'Coding','Non-coding'},'Location','best');
    
    grid on


    %% ==========================================
    %% TEMPERATURE
    %% ==========================================
    
    nexttile
    
    semilogy(result.history.temperature,'m','LineWidth',1.5);
    
    xlabel('Iteration');
    ylabel('Temperature');
    
    title('Cooling schedule');
    
    grid on


    %% ==========================================
    %% PURITY OF SUBPOP
    %% ==========================================

    nexttile
    
    purity = result.history.nCoding ./ result.history.size;

    plot(100*purity,'k','LineWidth',2);

    hold on
    
    xlabel('Iteration');
    
    ylabel('Purity');
    
    title('Population composition');

    grid on;


    %% ==========================================
    %% GLOBAL TITLE
    %% ==========================================
    
    sgtitle(sprintf('Simulated Annealing | Best P = %.3f | Best Population = [%s]',result.bestP,num2str(sort(result.bestPopulation))));

end