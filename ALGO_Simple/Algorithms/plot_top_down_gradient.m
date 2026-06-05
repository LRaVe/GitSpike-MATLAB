%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function plot_top_down_gradient(result,codingNeurons)

    N = length(result.populations{1});

    %% ==========================================
    %% Heatmap
    %% ==========================================

    M = nan(N,N);
    M(N,1) = result.P(1);

    for step = 1:length(result.candidateP)

        popSize = N-step;

        M(popSize,:) = result.candidateP{step};

    end

    %% ==========================================
    %% Figure
    %% ==========================================

    figure('Color','w');

    tiledlayout(1,2,'TileSpacing','compact');

    %% ==========================================
    %% LEFT PANEL
    %% ==========================================

    nexttile

    imagesc(M)

    set(gca,'YDir','normal')

    ax = gca;

    drawnow

    oldUnits = ax.Units;
    ax.Units = 'pixels';
    
    posPix = ax.Position;
    
    ax.Units = oldUnits;
    
    cellPix = min(posPix(3),posPix(4))/N;

    markerSize = max(4,0.8*cellPix);
    lineWidth  = max(1,0.15*cellPix);

    colormap(jet)

    colorbar

    hold on

    title('Top-Down','FontWeight','bold')

    xlabel('Neuron Index')
    ylabel('Size of population')

    xticks(1:N)
    yticks(1:N)

    %% --------------------------
    %% Deleted neurons
    %% --------------------------

    for step = 1:length(result.removedNeuron)

        popSize = N-step;

        neuron = result.removedNeuron(step);

        plot(neuron,popSize,'k_','MarkerSize',markerSize,'LineWidth',lineWidth)

    end

    %% --------------------------
    %% Complete population
    %% --------------------------

    plot(1,N,'ko','MarkerFaceColor','k','MarkerSize',0.5*markerSize,'LineWidth',lineWidth)

    %% --------------------------
    %% Best population
    %% --------------------------

    bestSize = length(result.bestPopulation);

    rectangle('Position',[0.5 bestSize-0.5 length(codingNeurons) 1],'EdgeColor',[0 1 0],'LineWidth',lineWidth)

    %% magenta cross

    for n = codingNeurons

        plot(n,bestSize,'mx','MarkerSize',1.2*markerSize,'LineWidth',lineWidth)

    end

    %% separation line between coding and non-coding

    xline(length(codingNeurons)+0.5,'k','LineWidth',4);


    %% ==========================================
    %% RIGHT PANEL
    %% ==========================================

    nexttile

    popSizes = cellfun(@length,result.populations);

    plot(result.P,popSizes,'-ok','LineWidth',2,'MarkerFaceColor','k')

    hold on

    plot(result.bestP,bestSize,'mo','MarkerSize',14,'LineWidth',3)

    plot(result.bestP,bestSize,'mx','MarkerSize',18,'LineWidth',3)

    xlabel('Best performance per pop size')
    ylabel('')

    ylim([0.5 N+0.5])

    grid on

end

