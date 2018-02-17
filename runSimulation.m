% in silico evolution of gene regulatory networks
% based on algorithm presented in "Design of genetic networks with
% specified functions by evolution in silico" by Francois and Hakim, PNAS
% (2004)
% Code written August 2017 as part of the MBL Physical Biology course
% Karna Gowda, Haneul Yoo, Ginger Hunter, Madhav Mani

nCells = 64;
nCycles = 400;
selPer = 10;

%initialize
for j = 1:nCells/2
    cycle(1).cell(j) = initialize(2); %initialize with two proteins
    cycle(1).cell(j+nCells/2) = cycle(1).cell(j);
end

%calculate initial fitness of first nCells/2 networks
for j = 1:nCells/2
    %simulate to obtain steady state
    [t,A] = ode15s(@(t,A) cycle(1).cell(j).K0(:) + cycle(1).cell(j).K1*A +...
        (cycle(1).cell(j).K2*A).*A + NL(A,cycle(1).cell(j).K3) + ...
        MM(A,cycle(1).cell(j).K4), [0 50], cycle(1).cell(j).A0);
    %perturb a little and run on a uniform grid
    [t,A] = ode15s(@(t,A) cycle(1).cell(j).K0(:) + cycle(1).cell(j).K1*A +...
        (cycle(1).cell(j).K2*A).*A + NL(A,cycle(1).cell(j).K3) +...
        MM(A,cycle(1).cell(j).K4), 0:0.5:100, [1.05*A(end,1) A(end,2:end)]);

    fit(1,j) = fitness(A(:,1),t(end),selPer,false); %calculate fitness
    netsize(1,j) = length(cycle(1).cell(j).A0); %compute number of variables in network
end

f=figure('position',[100 100 1024 640]);

%mutate & select
fitmean = nan(nCycles,1);
ancestor = nan(nCycles,nCells);
fitPlot = [];
for i = 1:nCycles
    tic   
    for j = nCells/2+1:nCells %mutate and calculate fitness for the second half of cells
        cycle(i).cell(j) = mutate(cycle(i).cell(j), 2); %apply two mutations to each cell
        %simulate to obtain steady state
        [t,A] = ode15s(@(t,A) cycle(i).cell(j).K0(:) + cycle(i).cell(j).K1*A +...
            (cycle(i).cell(j).K2*A).*A + NL(A,cycle(i).cell(j).K3) + ...
            MM(A,cycle(i).cell(j).K4), [0 50], cycle(i).cell(j).A0);
%         %perturb a little and run on a uniform grid
        [t,A] = ode15s(@(t,A) cycle(i).cell(j).K0(:) + cycle(i).cell(j).K1*A +...
            (cycle(i).cell(j).K2*A).*A + NL(A,cycle(i).cell(j).K3) +...
            MM(A,cycle(i).cell(j).K4), 0:0.5:100, [1.05*A(end,1) A(end,2:end)]);
        
        fit(i,j) = fitness(A(:,1),t(end),selPer,false); %calculate fitness
        netsize(i,j) = length(cycle(i).cell(j).A0);

        if sum(A(end,:)<-1e-4) > 0
            disp('negative solution detected')
            fit(i,j) = NaN;
        end
        if sum(A(end,:)>1e4) > 0
            disp('diverging solution detected')
            fit(i,j) = NaN;
        end
    end
    
    netmean(i) = mean(netsize(i,:)); %compute mean network size
    fitTemp = fit(i,:); 
    fitmean(i) = nanmean(fitTemp); %compute mean fitness over all cells
    fitTemp(isnan(fitTemp)) = -Inf; %send NaN values of fitness to -Inf for sorting
    [s,idx] = sort(fitTemp,'descend'); %sort cells by fitness
    
    %plotting
    subplot(2,21,1:5)
    %histogram(s,10)
    hist(s,10)
    xlabel('fitness')
    ylabel('frequency')
    ylim([0 nCells])
    subplot(2,21,22:30)
    plot(1:nCycles,fitmean,'r','LineWidth',2)
    xlim([0 nCycles])
    ylim([0 1])
    xlabel('cycle number')
    ylabel('mean fitness')
    windows = [8:11;13:16;18:21];
    plotCells = [1 nCells/2 find(isfinite(s) ,1,'last')];
    for k = 1:3
        subplot(2,21,windows(k,:))
        jj = idx(plotCells(k));
        %simulate to obtain steady state
        [t,A] = ode15s(@(t,A) cycle(i).cell(jj).K0(:) + cycle(i).cell(jj).K1*A +...
            (cycle(i).cell(jj).K2*A).*A + NL(A,cycle(i).cell(jj).K3) +...
            MM(A,cycle(i).cell(jj).K4), [0 100], cycle(i).cell(jj).A0);
        %perturb a little and run on a uniform grid
        [t,A] = ode15s(@(t,A) cycle(i).cell(jj).K0(:) + cycle(i).cell(jj).K1*A +...
            (cycle(i).cell(jj).K2*A).*A + NL(A,cycle(i).cell(jj).K3) +...
            MM(A,cycle(i).cell(jj).K4), 0:0.5:100, [1.05*A(end,1) A(end,2:end)]);
        
        plot(t,A(:,1),'k','LineWidth',2)
        if k == 1
        ylabel('[A]')
        end
        xlabel('time')
        title(['fitness = ', num2str(s(plotCells(k))),', N = ',num2str(size(A,2))])
    end
    subplot(2,21,33:42)
    fitPlot = [fitPlot; s];
    imagesc(fitPlot')
    xlim([0 nCycles])
    set(gca,'YDir','reverse')
    xlabel('cycle')
    ylabel('fitness rank')
    h = colorbar;
    ylabel(h, 'fitness')
    caxis([0 1])
    drawnow
    
    %propagate fittest cells to the next round
    cellTemp(1) = cycle(1).cell(1);
    cellTemp(nCells) = cycle(1).cell(1);
    for j = 1:nCells/2
        cellTemp(j) = cycle(i).cell(idx(j));
        ancestor(i,j) = idx(j); %store index of ancestor from this cycle
        cellTemp(j+nCells/2) = cycle(i).cell(idx(j));
        ancestor(i,j+nCells/2) = idx(j); %store index of ancestor from this cycle
        fit(i+1,j) = s(j);
        netsize(i+1,j) = netsize(i,idx(j));
        if i == 1
            ID{1,j} = ['0_',sprintf('%05d',j)];
            ID{1,j+nCells/2} = ['1_',sprintf('%05d',j)];
        else
            ID{i,j}             = ['0',ID{i-1,idx(j)}];
            ID{i,j+nCells/2}    = ['1',ID{i-1,idx(j)}];
        end
    end
    cycle(i+1).cell = cellTemp;
    runtime(i) = toc; %store run time for this cycle
    runtime(i)
end