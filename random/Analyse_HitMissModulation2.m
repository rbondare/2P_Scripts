%%

climStart = 0;
climStop = 10;
cmap = EntropyColourmap;
DotSize = 80;
conditions = {'Hit', 'Miss'};
fig1 = figure('Position', [134 160 1311 675]);
tt = tiledlayout(6,4);
for trial = 1:12
    % Plot mean activity
    nexttile;
    Periods = [83:85];
    [c, I] = sort(mean(Activity.Position1(:, Periods, Performance.hit(trial)), [2,3]));
    sc = scatter(CaData.Ca_centroid_voxel(I, 4), CaData.Ca_centroid_voxel(I, 5), DotSize, c, 'filled', 'MarkerEdgeColor', [0.6 0.6 0.6]);
    sc.LineWidth = 0.2;
    axis equal; axis off;
    clim([climStart climStop]);   
        colormap(cmap);
        if trial < 3
            subtitle('Onset')
        end
         nexttile;

             Periods = [83:105];
          [~,a] = max(mean(Activity.Position1(:, Periods,  Performance.hit(trial),1)),[],2);
    [c, I] = sort(mean(Activity.Position1(:, 82+a,  Performance.hit(trial)), [2,3]));
    sc = scatter(CaData.Ca_centroid_voxel(I, 4), CaData.Ca_centroid_voxel(I, 5), DotSize, c, 'filled', 'MarkerEdgeColor', [0.6 0.6 0.6]);
    sc.LineWidth = 0.2;
    axis equal; axis off;
    clim([climStart climStop]);   
        colormap(cmap);
             if trial < 3
            subtitle('Strongest Frame')
        end
end
% Add a single colorbar for the entire figure
colorbar('Position', [0.92 0.11 0.02 0.815]);
% obj = scalebar;obj.XUnit = '\mum';obj.YLen = 0;obj.hTextY.String  =[];obj.Position = [55, 30];


%%
fig1 = figure('Position', [134 160 1311 675]);
tt = tiledlayout(6,4);
for trial = 1:12
    % Plot mean activity
    nexttile;
    Periods = [83:85];
    [c, I] = sort(mean(Activity.Position1(:, Periods, Performance.miss(trial)), [2,3]));
    sc = scatter(CaData.Ca_centroid_voxel(I, 4), CaData.Ca_centroid_voxel(I, 5), DotSize, c, 'filled', 'MarkerEdgeColor', [0.6 0.6 0.6]);
    sc.LineWidth = 0.2;
    axis equal; axis off;
    clim([climStart climStop]);   
        colormap(cmap);
        if trial < 3
            subtitle('Onset')
        end
         nexttile;

             Periods = [83:105];
          [~,a] = max(mean(Activity.Position1(:, Periods,  Performance.miss(trial),1)),[],2);
    [c, I] = sort(mean(Activity.Position1(:, 82+a,  Performance.miss(trial)), [2,3]));
    sc = scatter(CaData.Ca_centroid_voxel(I, 4), CaData.Ca_centroid_voxel(I, 5), DotSize, c, 'filled', 'MarkerEdgeColor', [0.6 0.6 0.6]);
    sc.LineWidth = 0.2;
    axis equal; axis off;
    clim([climStart climStop]);   
        colormap(cmap);
             if trial < 3
            subtitle('Strongest Frame')
        end
end
% Add a single colorbar for the entire figure
colorbar('Position', [0.92 0.11 0.02 0.815]);
% obj = scalebar;obj.XUnit = '\mum';obj.YLen = 0;obj.hTextY.String  =[];obj.Position = [55, 30];


%%
responseNumberHit = [];
responseNumberMiss = [];
for r = 1:size(Rec,1)
    tempHit = [];
    tempMiss = [];
    for t = Stimuli(r).TrialsPosition1
        NotGoodSize =1;
        Counter = 0;
        while NotGoodSize & Counter < 3
        [idx, C] = kmeans(Activity(r).allZ(:,82:130,t),2);
        if  length(find(idx == 1)) > 50 & length(find(idx == 2)) > 40
            NotGoodSize = 0;
        end
                    Counter = Counter + 1;
        end
        if length(find(idx == 1)) > length(find(idx == 2))
            ClusterID = 2;
        else
            ClusterID = 1;
        end
        if ismember(t,Performance(r).hit)
            tempHit = [tempHit, length(find(idx == ClusterID))];
        else
            tempMiss = [tempMiss, length(find(idx == ClusterID))];
        end

    end
    responseNumberHit = [responseNumberHit,{tempHit}];
    responseNumberMiss = [responseNumberMiss,{tempMiss}];
    % [mean(responseNumberHit),mean(responseNumberMiss); length(responseNumberHit) ,length(responseNumberMiss)]
end
responseNumberHit
responseNumberHit_All = horzcat(responseNumberHit{:})';
responseNumberMiss_All = horzcat(responseNumberMiss{:})';
disp(join(["Neuron number responding to stimulus during hit: ", string((mean(responseNumberHit_All)))]))
disp(join(["Neuron number responding to stimulus during miss: ", string((mean(responseNumberMiss_All)))]))
data = [responseNumberHit_All;responseNumberMiss_All];
xl = cellstr([repmat("Hit",length(responseNumberHit_All),1);repmat("Miss",length(responseNumberMiss_All),1)]);
xl = categorical(xl);
fig=figure('Position',[ 134         160        1311         675]);
vs = violinplot(data, xl ,'ShowData',true, 'ViolinAlpha', 0.5, 'ShowMean', true, 'EdgeColor' , [0.1 0.1 0.1],  'ViolinColor', [0.1 0.1 0.1]);
xlabel("Performance", 'FontSize',16,'FontWeight','bold')
ylabel('Size of stimulus Cluster', 'FontSize',16,'FontWeight','bold')
fontsize(24,"points")

fig=figure('Position',[ 134         160        1311         675]);
f = tiledlayout(1,length(responseNumberHit));
for a = 1:length(responseNumberHit)
    ax(a) = nexttile;
    data = [responseNumberHit{1,a}';responseNumberMiss{1,a}'];
    xl = cellstr([repmat("Hit",length(responseNumberHit{1,a}),1);repmat("Miss",length(responseNumberMiss{1,a}),1)]);
    xl = categorical(xl);
    vs = violinplot(data, xl ,'ShowData',true, 'ViolinAlpha', 0.5, 'ShowMean', true, 'EdgeColor' , [0.1 0.1 0.1],  'ViolinColor', [0.1 0.1 0.1]);
    xlabel("Performance", 'FontSize',16,'FontWeight','bold')
    ylabel('Size of stimulus Cluster', 'FontSize',16,'FontWeight','bold')
end
linkaxes([ax],'xy')
fontsize(24,"points")

if ~isempty(HitControl)
responseNumberHitControl = [];
responseNumberMissControl = [];
for r = 1:size(RecControl,1)
    tempHit = [];
    tempMiss = [];
    for t = StimuliControl(r).TrialsPosition1

         NotGoodSize =1;
        Counter = 0;
        while NotGoodSize & Counter < 3
        [idx, C] = kmeans(ActivityControl(r).allZ(:,82:130,t),2);
        if  length(find(idx == 1)) > 50 & length(find(idx == 2)) > 40
            NotGoodSize = 0;
        end
                    Counter = Counter + 1;
        end
        if length(find(idx == 1)) > length(find(idx == 2))
            ClusterID = 2;
        else
            ClusterID = 1;
        end
        if ismember(t,PerformanceControl(r).hit)
            tempHit = [tempHit, length(find(idx == ClusterID))];
        else
            tempMiss = [tempMiss, length(find(idx == ClusterID))];
        end

    end
    responseNumberHitControl = [responseNumberHitControl,{tempHit}];
    responseNumberMissControl = [responseNumberMissControl,{tempMiss}];
    % [mean(responseNumberHit),mean(responseNumberMiss); length(responseNumberHit) ,length(responseNumberMiss)]
end
responseNumberHit_AllControl = horzcat(responseNumberHitControl{:})';
responseNumberMiss_AllControl = horzcat(responseNumberMissControl{:})';
disp(join(["Neuron number responding to stimulus during hit: ", string((mean(responseNumberHit_AllControl)))]))
disp(join(["Neuron number responding to stimulus during miss: ", string((mean(responseNumberMiss_AllControl)))]))

data = [responseNumberHit_AllControl;responseNumberMiss_AllControl];
xl = cellstr([repmat("Hit",length(responseNumberHit_AllControl),1);repmat("Miss",length(responseNumberMiss_AllControl),1)]);
xl = categorical(xl);
fig=figure('Position',[ 134         160        1311         675]);
vs = violinplot(data, xl ,'ShowData',true, 'ViolinAlpha', 0.5, 'ShowMean', true, 'EdgeColor' , [0.1 0.1 0.1],  'ViolinColor', [0.1 0.1 0.1]);
xlabel("Performance", 'FontSize',16,'FontWeight','bold')
ylabel('Size of stimulus Cluster', 'FontSize',16,'FontWeight','bold')
fontsize(24,"points")


fig=figure('Position',[ 134         160        1311         675]);
f = tiledlayout(1,length(responseNumberHitControl));
for a = 1:length(responseNumberHitControl)
    ax(a) = nexttile;
    data = [responseNumberHitControl{1,a}';responseNumberMissControl{1,a}'];
    xl = cellstr([repmat("Hit",length(responseNumberHitControl{1,a}),1);repmat("Miss",length(responseNumberMissControl{1,a}),1)]);
    xl = categorical(xl);
    vs = violinplot(data, xl ,'ShowData',true, 'ViolinAlpha', 0.5, 'ShowMean', true, 'EdgeColor' , [0.1 0.1 0.1],  'ViolinColor', [0.1 0.1 0.1]);
    xlabel("Performance", 'FontSize',16,'FontWeight','bold')
    ylabel('Size of stimulus Cluster', 'FontSize',16,'FontWeight','bold')
end
linkaxes([ax],'xy')
fontsize(24,"points")
end

if ~isempty(HitSnapshot)
responseNumberHitSnapshot = [];
responseNumberMissSnapshot = [];
for r = 1:size(RecSnapshot,1)
    tempHit = [];
    tempMiss = [];
    for t = StimuliSnapshot(r).TrialsPosition1

NotGoodSize =1;
        Counter = 0;
        while NotGoodSize & Counter < 3
        [idx, C] = kmeans(ActivitySnapshot(r).allZ(:,52:80,t),2);
        if  length(find(idx == 1)) > 50 & length(find(idx == 2)) > 40
            NotGoodSize = 0;
        end
                    Counter = Counter + 1;
        end


        if length(find(idx == 1)) > length(find(idx == 2))
            ClusterID = 2;
        else
            ClusterID = 1;
        end
        % if ismember(t,PerformanceSnapshot(r).hit)
        %     tempHit = [tempHit, length(find(idx == ClusterID))];
        % else
            tempMiss = [tempMiss, length(find(idx == ClusterID))];
        % end

    end
    responseNumberHitSnapshot = [responseNumberHitSnapshot,{tempHit}];
    responseNumberMissSnapshot = [responseNumberMissSnapshot,{tempMiss}];
    % [mean(responseNumberHit),mean(responseNumberMiss); length(responseNumberHit) ,length(responseNumberMiss)]
end
responseNumberHit_AllSnapshot = horzcat(responseNumberHitSnapshot{:})';
responseNumberMiss_AllSnapshot = horzcat(responseNumberMissSnapshot{:})';
disp(join(["Neuron number responding to stimulus during hit: ", string((mean(responseNumberHit_AllSnapshot)))]))
disp(join(["Neuron number responding to stimulus during miss: ", string((mean(responseNumberMiss_AllSnapshot)))]))

data = [responseNumberHit_AllSnapshot;responseNumberMiss_AllSnapshot];
xl = cellstr([repmat("Hit",length(responseNumberHit_AllSnapshot),1);repmat("Miss",length(responseNumberMiss_AllSnapshot),1)]);
xl = categorical(xl);
fig=figure('Position',[ 134         160        1311         675]);
vs = violinplot(data, xl ,'ShowData',true, 'ViolinAlpha', 0.5, 'ShowMean', true, 'EdgeColor' , [0.1 0.1 0.1],  'ViolinColor', [0.1 0.1 0.1]);
xlabel("Performance", 'FontSize',16,'FontWeight','bold')
ylabel('Size of stimulus Cluster', 'FontSize',16,'FontWeight','bold')

fig=figure('Position',[ 134         160        1311         675]);
f = tiledlayout(1,length(responseNumberHitSnapshot));
for a = 1:length(responseNumberHitSnapshot)
    ax(a) = nexttile;
    data = [responseNumberHitSnapshot{1,a}';responseNumberMissSnapshot{1,a}'];
    xl = cellstr([repmat("Hit",length(responseNumberHitSnapshot{1,a}),1);repmat("Miss",length(responseNumberMissSnapshot{1,a}),1)]);
    xl = categorical(xl);
    vs = violinplot(data, xl ,'ShowData',true, 'ViolinAlpha', 0.5, 'ShowMean', true, 'EdgeColor' , [0.1 0.1 0.1],  'ViolinColor', [0.1 0.1 0.1]);
    xlabel("Performance", 'FontSize',16,'FontWeight','bold')
    ylabel('Size of stimulus Cluster', 'FontSize',16,'FontWeight','bold')
end
linkaxes([ax],'xy')
fontsize(24,"points")
end
%%
data = [responseNumberHit_All;responseNumberMiss_All;rot90(responseNumberMissSnapshot{:})];
xl = cellstr([repmat("Hit",length(responseNumberHit_All),1);repmat("Miss",length(responseNumberMiss_All),1);repmat("NoSpout",length(responseNumberMissSnapshot{:}),1)]);
xl = categorical(xl);
fig=figure('Position',[ 134         160        1311         675]);
vs = violinplot(data, xl ,'ShowData',true, 'ViolinAlpha', 0.5, 'ShowMean', true, 'EdgeColor' , [0.1 0.1 0.1],  'ViolinColor', [0.1 0.1 0.1]);
xlabel("Performance", 'FontSize',16,'FontWeight','bold')
ylabel('Size of stimulus Cluster', 'FontSize',16,'FontWeight','bold')
fontsize(24,"points")

%         data = {};
%         data{1,1} =responseNumberHit';
%             data{1,2} = responseNumberMiss';
%
%    names = ['Expert']
%     c =[params.colour(6,:);params.colour(7,:)];
% group_inx = [repmat([repmat(1,size(data{1},1),1) ; repmat(2,size(data{2},1),1)],1,1) ];
%     group_names = {'Hit','Miss'};
%
%     fig=figure('Position',[ 134         160        1311         675]);
% % half-violin plots with white boxplots, jittered scatter and linkline
% h = daviolinplot(data,'groups',group_inx,'colors',c,'box',3,...
%     'boxcolor','same','scatter',2,'jitter',1,'scattercolor','same',...
%     'scattersize',10,'scatteralpha',1,'linkline',1, 'withinlines',0,...
%     'xtlabels', names,'violin','half', 'legend',group_names);
%
% ylabel('Stimulus cluster size');
% xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]); % make more space for the legend
% set(h.sc,'MarkerEdgeColor','none');      % remove marker edge color
% set(gca,'FontSize',10);    hold on
%     legend boxoff



%%
climStart = -3;
climStop = 3;

% climStart = 1;
% climStop = 4;

showFigures = 0;
DotSize = 400;
MomentOfStrongestActivity = [];

for r = 1:size(Rec,1)
    [M,I] = max(mean(Hit.AllTrials{r,1}(:,Time.StimStart+3:Time.StimStop,:),[1,3])./mean(Miss.AllTrials{r,1}(:,Time.StimStart+3:Time.StimStop,:),[1,3]));
    for Trial = 1:size(Hit.AllTrials{r,1},3)
        [M,Itrial] =  max(Hit.AllTrials{r,1}(:,80:100,Trial),[],2);
        MomentOfStrongestActivity = [MomentOfStrongestActivity, Itrial/10];
        figure;histogram(Itrial/10,'BinWidth',0.1); xline(Hit.ReactionTimes(Trial),'LineWidth',3,'Color','r')
    end
        % figure;title(Itrial);
        % [c,Isort] = sort(mean(Hit.AllTrials{r,1}(:,80+Itrial,Trial),2));
        % sc = scatter(CaData(r).Ca_centroid_voxel(Isort,4),CaData(r).Ca_centroid_voxel(Isort,5),DotSize,c,'filled','MarkerEdgeColor',[0.6 0.6 0.6]);
        % sc(1).LineWidth = 0.2;
        % axis equal;
        % axis off;
        % clim([climStart climStop]);
        if showFigures

            for Trial = 1:size(Hit.AllTrials{r,1},3)
        [M,Itrial] =  max(Hit.AllTrials{r,1}(:,80:120,Trial),[],2);
        fig=figure('Position',[ 134         160        1311         675]);     
        [c,Isort] = sort(mean(Hit.AllCells_AllTrials{r,1}(:,80+Itrial,Trial),2));
        sc = scatter(CaData(r).Ca_centroid_voxel(Isort,4),CaData(r).Ca_centroid_voxel(Isort,5),DotSize,c,'filled','MarkerEdgeColor',[0.6 0.6 0.6]);
       title('Hit');
        sc(1).LineWidth = 0.2;
        axis equal;
        axis off;
       clim([climStart climStop]);
       colormap(params.cmap)

    end

     for Trial = 1:size(Miss.AllTrials{r,1},3)
               
    fig=figure('Position',[ 134         160        1311         675]);     
        [c,Isort] = sort(mean(Miss.AllCells_AllTrials{r,1}(:,80+Itrial,Trial),2));
        sc = scatter(CaData(r).Ca_centroid_voxel(Isort,4),CaData(r).Ca_centroid_voxel(Isort,5),DotSize,c,'filled','MarkerEdgeColor',[0.6 0.6 0.6]);
       title('Miss');
        sc(1).LineWidth = 0.2;
        axis equal;
        axis off;
       clim([climStart climStop]);
       colormap(params.cmap)
     end
        end



     climStart = 1;
climStop = 4;
    fig=figure('Position',[ 134         160        1311         675]);
    f = tiledlayout(ceil(sqrt(size(Hit.AllTrials{r,1},3))),ceil(sqrt(size(Hit.AllTrials{r,1},3))));
    for Trial = 1:size(Hit.AllTrials{r,1},3)
        [M,FirstTrial] = max(mean(Hit.AllTrials{r,1}(:,Time.StimStart+3:Time.StimStop,Trial),[1]));
    [M,CellOrder] =  max(Hit.AllTrials{r,1}(:,FirstTrial+3,1),[],2);
        nexttile;
        [M,Itrial] =  max(Hit.AllTrials{r,1}(CellOrder,80:100,Trial),[],2);
        imagesc(Hit.AllTrials{r,1}(:,70:120,Trial))
        % clim([climStart climStop]);
        colormap(params.cmap)
    end

    fig=figure('Position',[ 134         160        1311         675]);
    f = tiledlayout(ceil(sqrt(size(Hit.AllTrials{r,1},3))),ceil(sqrt(size(Hit.AllTrials{r,1},3))));
    [M,FirstTrial] = max(mean(Hit.AllTrials{r,1}(:,Time.StimStart+3:Time.StimStop,Trial)./Hit.AllTrials{r,1}(:,82,Trial),1));
    for Trial = 1:size(Hit.AllTrials{r,1},3)
        nexttile;
        [M,Itrial] =  max(Hit.AllTrials{r,1}(CellOrder,80:100,Trial),[],2);
        imagesc(Hit.AllTrials{r,1}(:,70:120,Trial)./Hit.AllTrials{r,1}(:,82,Trial))
         clim([climStart climStop]);
        colormap(params.cmap)
    end
 fig=figure('Position',[ 134         160        1311         675]);
    f = tiledlayout(3,3);
    [M,FirstTrial] = max(mean(Hit.AllTrials{r,1}(:,Time.StimStart+3:Time.StimStop,Trial)./Hit.AllTrials{r,1}(:,82,Trial),1));
    for Trial = 1:9
        nexttile;
        [M,Itrial] =  max(Hit.AllTrials{r,1}(CellOrder,80:100,Trial),[],2);
        imagesc(Hit.AllTrials{r,1}(:,70:120,Trial)./Hit.AllTrials{r,1}(:,82,Trial))
           clim([-6 6]);
        colormap(params.cmap)
         xlabel("Seconds after Stimulus")
         ylabel('Neuron ID')
    end
    fontsize(24,'points')


    fig=figure('Position',[ 134         160        1311         675]);
    f = tiledlayout(ceil(sqrt(size(Hit.AllTrials{r,1},3))),ceil(sqrt(size(Hit.AllTrials{r,1},3))));
    [M,FirstTrial] = max(mean(Hit.AllTrials{r,1}(:,Time.StimStart+3:Time.StimStop,Trial)./Hit.AllTrials{r,1}(:,82,Trial),1));
    [M,CellOrder] =  sort(Hit.AllTrials{r,1}(:,FirstTrial+3,1)./Hit.AllTrials{r,1}(:,82,Trial),1);
    for Trial = 1:size(Hit.AllTrials{r,1},3)
        nexttile;
        imagesc(Hit.AllTrials{r,1}(CellOrder,70:120,Trial)./Hit.AllTrials{r,1}(CellOrder,82,Trial))
         clim([-1 climStop]);
        colormap(params.cmap)
    end

     fig=figure('Position',[ 134         160        1311         675]);
    f = tiledlayout(3,3);
    [M,FirstTrial] = max(mean(Hit.AllTrials{r,1}(:,Time.StimStart+3:Time.StimStop,Trial)./Hit.AllTrials{r,1}(:,82,Trial),1));
    [M,CellOrder] =  sort(Hit.AllTrials{r,1}(:,FirstTrial+3,1)./Hit.AllTrials{r,1}(:,82,Trial),1);
    for Trial = 1:9
        nexttile;
        imagesc(linspace(0,4,40),1:length(CellOrder),Hit.AllTrials{r,1}(CellOrder,80:120,Trial)./Hit.AllTrials{r,1}(CellOrder,82,Trial))
          clim([-6 6]);
        colormap(params.cmap)
         xlabel("Seconds after Stimulus")
         ylabel('Neuron ID')
    end
    fontsize(24,'points')
   

end

%%

%%

low= [];high = [];L = 0;
for r = 1:size(Rec,1)
    input = mean(Activity(r).all(clusters(r).StimCluster,Timing(r).StimStartFrameNum:Timing(r).StimEndFrameNum,Stimuli{1, r}.TrialsPosition1),3);
    halfMax = max(input,[],'all')/2;
    temp_high =  find ((max(input,[],[2,3]) > halfMax) == 1);
    temp_low = find ((max(input,[],[2,3]) > halfMax) == 0);
    high= [high;temp_high + L]; low = [low; temp_low + L]; L = L + size(input, 1);
end

figure;scatter(mean(Hit.Stim_P1(low,Time.StimStart+I-1),3),mean(Miss.Stim_P1(low,Time.StimStart+I-1),3),[],params.colour(2,:),'filled');hold on
scatter(mean(Hit.Stim_P1(high,Time.StimStart+I-1),3),mean(Miss.Stim_P1(high,Time.StimStart+I-1),3),[],params.colour(5,:),'filled');
hline1 = refline(1); hline1.LineWidth = 3 ;  hline1.Color = 'r'; ax1 = gca;  hline =  lsline;
hline(1).LineWidth = 2; hline(1).Color = params.colour(5,:);
hline(2).LineWidth = 2;  hline(2).Color = params.colour(2,:);
xlabel(' delta F/F hit trials');
ylabel('delta F/F miss trials');
title('Clustered based on maximum activity in response to stimulus')

low= [];med = []; high = [];L = 0;
for r = 1:size(Rec,1)
    input = mean(Activity(r).all(clusters(r).StimCluster,Timing(r).StimStartFrameNum:Timing(r).StimEndFrameNum,Stimuli{1, r}.TrialsPosition1),3);
    twothirdMax = max(input,[],'all')/3*2;
    onethirdMax = max(input,[],'all')/3;
    temp_high =  find ((max(input,[],[2,3]) > twothirdMax) == 1);
    temp_med = find ((max(input,[],[2,3]) <  twothirdMax) == 1 & (max(input,[],[2,3]) >  onethirdMax) == 1 );
    temp_low = find ((max(input,[],[2,3]) < onethirdMax) == 1);
    high= [high;temp_high + L]; med = [med; temp_med + L ]; low = [low; temp_low + L]; L = L + size(input, 1);
end

figure;scatter(mean(Hit.Stim_P1(low,Time.StimStart+I-1),3),mean(Miss.Stim_P1(low,Time.StimStart+I-1),3),[],params.colour(2,:),'filled');hold on
scatter(mean(Hit.Stim_P1(med,Time.StimStart+I-1),3),mean(Miss.Stim_P1(med,Time.StimStart+I-1),3),[],params.colour(4,:),'filled');
scatter(mean(Hit.Stim_P1(high,Time.StimStart+I-1),3),mean(Miss.Stim_P1(high,Time.StimStart+I-1),3),[],params.colour(5,:),'filled');
hline1 = refline(1); hline1.LineWidth = 3 ;  hline1.Color = 'r'; ax1 = gca;  hline =  lsline;
hline(1).LineWidth = 2; hline(1).Color = params.colour(5,:);
hline(2).LineWidth = 2; hline(2).Color = params.colour(4,:);
hline(3).LineWidth = 2;  hline(3).Color = params.colour(2,:);
xlabel(' delta F/F hit trials');
ylabel('delta F/F miss trials');
title('Clustered based on above/belove max response to stimulus')



low= [];; high = [];L = 0;
for r = 1:size(Rec,1)
    input = mean(Activity(r).all(clusters(r).StimCluster,Timing(r).StimStartFrameNum:Timing(r).StimEndFrameNum,Stimuli{1, r}.TrialsPosition1),[2,3]);
    tempMean = mean(input,'all');
    temp_high =  find ((mean(input,[2,3]) > tempMean) == 1);
    temp_low = find ((mean(input,[2,3]) < tempMean) == 1);
    high= [high;temp_high + L]; low = [low; temp_low + L]; L = L + size(input, 1);
end

figure;scatter(mean(Hit.Stim_P1(low,Time.StimStart+I-1),3),mean(Miss.Stim_P1(low,Time.StimStart+I-1),3),100,params.colour(2,:),'filled');hold on
scatter(mean(Hit.Stim_P1(high,Time.StimStart+I-1),3),mean(Miss.Stim_P1(high,Time.StimStart+I-1),3),100,params.colour(5,:),'filled');
hline1 = refline(1); hline1.LineWidth = 3 ;  hline1.Color = 'r'; ax1 = gca;
hline =  lsline;hline(1).LineWidth = 2; hline(1).Color = params.colour(5,:);
hline(2).LineWidth = 2;  hline(2).Color = params.colour(2,:);
xlabel(' delta F/F hit trials');
ylabel('delta F/F miss trials');
title('Clustered based on above/belove mean response to stimulus')

low= [];; high = [];L = 0;
for r = 1:size(Rec,1)
    input = mean(Activity(r).all(clusters(r).StimCluster,Timing(r).StimStartFrameNum:Timing(r).StimEndFrameNum,Stimuli{1, r}.TrialsPosition1),[2,3]);
    tempMedian = median(input,'all');
    temp_high =  find ((median(input,[2,3]) > tempMedian) == 1);
    temp_low = find ((median(input,[2,3]) < tempMedian) == 1);
    high= [high;temp_high + L]; low = [low; temp_low + L]; L = L + size(input, 1);
end

figure;scatter(mean(Hit.Stim_P1(low,Time.StimStart+I-1),3),mean(Miss.Stim_P1(low,Time.StimStart+I-1),3),100,params.colour(2,:),'filled');hold on
scatter(mean(Hit.Stim_P1(high,Time.StimStart+I-1),3),mean(Miss.Stim_P1(high,Time.StimStart+I-1),3),100,params.colour(5,:),'filled');
hline1 = refline(1); hline1.LineWidth = 3 ;  hline1.Color = 'r'; ax1 = gca;
hline =  lsline;hline(1).LineWidth = 2; hline(1).Color = params.colour(5,:);
hline(2).LineWidth = 2;  hline(2).Color = params.colour(2,:);
xlabel(' delta F/F hit trials');
ylabel('delta F/F miss trials');
title('Clustered based on above/belove median response to stimulus')

%%

nClusters = 2;
n = [];
for r = 1:size(Rec,1)
    input = Activity(r).all(clusters(r).StimCluster,Timing(r).StimStartFrameNum:Timing(r).StimEndFrameNum,Performance(r).hit);
    [tempn, C] = kmeans(reshape(input,size(input,1),size(input,2) * size(input,3)),nClusters);
    n = [n;tempn];
end
figure;scatter(mean(Hit.Stim_P1(find(n == 1),Time.StimStart+I-1),3),mean(Miss.Stim_P1(find(n == 1),Time.StimStart+I-1),3),100,params.colour(2,:),'filled');hold on
scatter(mean(Hit.Stim_P1(find(n == 2),Time.StimStart+I-1),3),mean(Miss.Stim_P1(find(n == 2),Time.StimStart+I-1),3),100,params.colour(5,:),'filled');
hline1 = refline(1); hline1.LineWidth = 3 ;  hline1.Color = 'r'; ax1 = gca;
hline =  lsline;hline(1).LineWidth = 2; hline(1).Color = params.colour(5,:);
hline(2).LineWidth = 2;  hline(2).Color = params.colour(2,:);
title('k=means clustered, hit trials, stimulus period')
xlabel(' delta F/F hit trials');
ylabel('delta F/F miss trials');


nClusters = 2;
n = [];
for r = 1:size(Rec,1)
    input = Activity(r).all(clusters(r).StimCluster,Timing(r).StimStartFrameNum:Timing(r).StimEndFrameNum,Stimuli{1,r}.TrialsPosition1);
    [tempn, C] = kmeans(reshape(input,size(input,1),size(input,2) * size(input,3)),nClusters);
    n = [n;tempn];
end
figure;scatter(mean(Hit.Stim_P1(find(n == 1),Time.StimStart+I-1),3),mean(Miss.Stim_P1(find(n == 1),Time.StimStart+I-1),3),100,params.colour(2,:),'filled');hold on
scatter(mean(Hit.Stim_P1(find(n == 2),Time.StimStart+I-1),3),mean(Miss.Stim_P1(find(n == 2),Time.StimStart+I-1),3),100,params.colour(5,:),'filled');
hline1 = refline(1); hline1.LineWidth = 3 ;  hline1.Color = 'r'; ax1 = gca;
hline =  lsline;hline(1).LineWidth = 2; hline(1).Color = params.colour(5,:);
hline(2).LineWidth = 2;  hline(2).Color = params.colour(2,:);
title('k=means clustered, all trials, stimulus period')
xlabel(' delta F/F hit trials');
ylabel('delta F/F miss trials');

nClusters = 2;
n = [];
for r = 1:size(Rec,1)
    input = Activity(r).all(clusters(r).StimCluster,:,Stimuli{1,r}.TrialsPosition1);
    [tempn, C] = kmeans(reshape(input,size(input,1),size(input,2) * size(input,3)),nClusters);
    n = [n;tempn];
end
figure;scatter(mean(Hit.Stim_P1(find(n == 1),Time.StimStart+I-1),3),mean(Miss.Stim_P1(find(n == 1),Time.StimStart+I-1),3),100,params.colour(2,:),'filled');hold on
scatter(mean(Hit.Stim_P1(find(n == 2),Time.StimStart+I-1),3),mean(Miss.Stim_P1(find(n == 2),Time.StimStart+I-1),3),100,params.colour(5,:),'filled');
hline1 = refline(1); hline1.LineWidth = 3 ;  hline1.Color = 'r'; ax1 = gca;
hline =  lsline;hline(1).LineWidth = 2; hline(1).Color = params.colour(5,:);
hline(2).LineWidth = 2;  hline(2).Color = params.colour(2,:);
title('k=means clustered, all trials, whole trial')
xlabel(' delta F/F hit trials');
ylabel('delta F/F miss trials');

%%
nClusters = 2;
n = [];
for r = 1:size(Rec,1)
    input = Activity(r).all(clusters(r).StimCluster,Timing(r).StimStartFrameNum:Timing(r).StimStartFrameNum+3,Stimuli{1,r}.TrialsPosition1);
    [tempn, C] = kmeans(reshape(input,size(input,1),size(input,2) * size(input,3)),nClusters);
    n = [n;tempn];
end

figure;scatter(mean(Hit.Stim_P1(find(n == 1),Time.StimStart+I-1),3),mean(Miss.Stim_P1(find(n == 1),Time.StimStart+I-1),3),100,params.colour(2,:),'filled');hold on
scatter(mean(Hit.Stim_P1(find(n == 2),Time.StimStart+I-1),3),mean(Miss.Stim_P1(find(n == 2),Time.StimStart+I-1),3),100,params.colour(5,:),'filled');
hline1 = refline(1); hline1.LineWidth = 3 ;  hline1.Color = 'r'; ax1 = gca;
hline =  lsline;hline(1).LineWidth = 2; hline(1).Color = params.colour(5,:);
hline(2).LineWidth = 2;  hline(2).Color = params.colour(2,:);
title('k=means clustered, all trials, stimulus onset')
xlabel(' delta F/F hit trials');
ylabel('delta F/F miss trials');
%%
nClusters = 2;
n = [];
for r = 1:size(Rec,1)
    input = Activity(r).all(clusters(r).StimCluster,Timing(r).StimEndFrameNum:Timing(r).StimEndFrameNum+2,Stimuli{1,r}.TrialsPosition1);
    [tempn, C] = kmeans(reshape(input,size(input,1),size(input,2) * size(input,3)),nClusters);
    n = [n;tempn];
end
figure;scatter(mean(Hit.Stim_P1(find(n == 1),Time.StimStart+I-1),3),mean(Miss.Stim_P1(find(n == 1),Time.StimStart+I-1),3),100,params.colour(2,:),'filled');hold on
scatter(mean(Hit.Stim_P1(find(n == 2),Time.StimStart+I-1),3),mean(Miss.Stim_P1(find(n == 2),Time.StimStart+I-1),3),100,params.colour(5,:),'filled');
hline1 = refline(1); hline1.LineWidth = 3 ;  hline1.Color = 'r'; ax1 = gca;
hline =  lsline;hline(1).LineWidth = 2; hline(1).Color = params.colour(5,:);
hline(2).LineWidth = 2;  hline(2).Color = params.colour(2,:);
title('k=means clustered, all trials, stimulus offset')
xlabel(' delta F/F hit trials');
ylabel('delta F/F miss trials');

%%
n = [];
for r = 1:size(Rec,1)
    on = mean(Activity(r).all(clusters(r).StimCluster,Timing(r).StimStartFrameNum:Timing(r).StimStartFrameNum+3,Stimuli{1,r}.TrialsPosition1),[2,3]);
    off = mean(Activity(r).all(clusters(r).StimCluster,Timing(r).StimEndFrameNum:Timing(r).StimEndFrameNum+3,Stimuli{1,r}.TrialsPosition1),[2,3]);
    tempn = ones(length(on),1);
    tempn(find( on - off > 0 )) = 2;
    n = [n;tempn];
end
figure;scatter(mean(Hit.Stim_P1(find(n == 1),Time.StimStart+I-1),3),mean(Miss.Stim_P1(find(n == 1),Time.StimStart+I-1),3),100,params.colour(2,:),'filled');hold on
scatter(mean(Hit.Stim_P1(find(n == 2),Time.StimStart+I-1),3),mean(Miss.Stim_P1(find(n == 2),Time.StimStart+I-1),3),100,params.colour(5,:),'filled');
hline1 = refline(1); hline1.LineWidth = 3 ;  hline1.Color = 'r'; ax1 = gca;
hline =  lsline;hline(1).LineWidth = 2; hline(1).Color = params.colour(5,:);
hline(2).LineWidth = 2;  hline(2).Color = params.colour(2,:);
title('On-response / Off-response')
xlabel(' delta F/F hit trials');
ylabel('delta F/F miss trials');



%%

figure;dscatter(mean(Hit.Stim_P1(:,Time.StimStart+I-1,:),3),mean(Miss.Stim_P1(:,Time.StimStart+I-1,:),3),'MSIZE',30)
hline1 = refline(1)
hline1.LineWidth = 5
% hold on
% hline = refline()
% hline.Color = 'r';
% hline.LineWidth = 5
xlabel(' delta F/F hit trials');
ylabel('delta F/F miss trials');

%%
HitRecording = [];
MissRecording = [];
C1 = 1;
C2 = length(clusters(1).StimCluster)  ;
for r = 1:size(Rec,1)
    HitRecording = [HitRecording;mean(Hit.Stim_P1(C1:C2,:),1)];
    MissRecording = [MissRecording;mean(Miss.Stim_P1(C1:C2,:),1)];
    C1 = C1 + length(clusters(r).StimCluster);
    C2 = C2 +  length(clusters(r).StimCluster)  ;
end

scatterX = HitRecording(:,Time.StimStart+I-1)
scatterY = MissRecording(:,Time.StimStart+I-1)
figure;dscatter(scatterX,scatterY ,'MSIZE',30)
hline1 = refline(1)
hline1.LineWidth = 5
% hold on
% hline = refline()
% hline.Color = 'r';
% hline.LineWidth = 5
xlabel(' delta F/F hit trials');
ylabel('delta F/F miss trials');

%%
HitSTD = [];
HitWithSTDHigh = [];
HitWithSTDLow = [];
AllWithSTDHigh = [];
AllWithSTDLow = [];
for r = 1:size(Rec,1)
    L =  length(HitSTD);
    tempSTD =std(Hit.AllTrials{r}(:,I,:),[],3);
    HitSTD = [HitSTD;  tempSTD];
    [~,Index] = sort(std(Hit.AllTrials{r}(:,I,:),[],3),1,'descend');
    HitWithSTDHigh =  [HitWithSTDHigh;Index(1:30)+L];
    HitWithSTDLow =  [HitWithSTDHigh;Index(31:end)+L];
    [~,Index] = sort(std(cat(3,Hit.AllTrials{r}(:,I,:),Miss.AllTrials{r}(:,I,:)),[],3),1,'descend');
    AllWithSTDHigh =  [AllWithSTDHigh;Index(1:30)+L];
    AllWithSTDLow =  [AllWithSTDHigh;Index(31:end)+L];

end
figure;hist(HitSTD,100)
title('distribution of standard deviation')


%%

nClusters = 2;
[n, C] = kmeans(HitSTD,nClusters);
figure;scatter(mean(Hit.Stim_P1(find(n == 1),Time.StimStart+I-1,:),3),mean(Miss.Stim_P1(find(n == 1),Time.StimStart+I-1,:),3));hold on
scatter(mean(Hit.Stim_P1(find(n == 2),Time.StimStart+I-1,:),3),mean(Miss.Stim_P1(find(n == 2),Time.StimStart+I-1,:),3));
hline1 = refline(1)
hline1.LineWidth = 3
hline1.Color = 'r';
hline =  lsline
hline(1).LineWidth = 2;
hline(1).Color = params.colour(1,:);
hline(2).LineWidth = 2;
hline(2).Color = params.colour(3,:);
xlabel(' delta F/F hit trials');
ylabel('delta F/F miss trials');


figure;scatter(mean(Hit.Stim_P1(HitWithSTDLow,Time.StimStart+I-1,:),3),mean(Miss.Stim_P1(HitWithSTDLow,Time.StimStart+I-1,:),3),100,params.colour(2,:),'filled');hold on
scatter(mean(Hit.Stim_P1(HitWithSTDHigh,Time.StimStart+I-1,:),3),mean(Miss.Stim_P1(HitWithSTDHigh,Time.StimStart+I-1,:),3),100,params.colour(5,:),'filled');
hline1 = refline(1); hline1.LineWidth = 3 ;  hline1.Color = 'r'; ax1 = gca;
hline =  lsline;hline(1).LineWidth = 2; hline(1).Color = params.colour(5,:);
hline(2).LineWidth = 2;  hline(2).Color = params.colour(2,:);
title('clustered based on high/low std, just hit trials')
legend('low','high')
xlabel(' delta F/F hit trials');
ylabel('delta F/F miss trials');

figure;scatter(mean(Hit.Stim_P1(AllWithSTDLow,Time.StimStart+I-1,:),3),mean(Miss.Stim_P1(AllWithSTDLow,Time.StimStart+I-1,:),3),100,params.colour(2,:),'filled');hold on
scatter(mean(Hit.Stim_P1(AllWithSTDHigh,Time.StimStart+I-1,:),3),mean(Miss.Stim_P1(AllWithSTDHigh,Time.StimStart+I-1,:),3),100,params.colour(5,:),'filled');
hline1 = refline(1); hline1.LineWidth = 3 ;  hline1.Color = 'r'; ax1 = gca;
hline =  lsline;hline(1).LineWidth = 2; hline(1).Color = params.colour(5,:);
hline(2).LineWidth = 2;  hline(2).Color = params.colour(2,:);
legend('low','high')
title('clustered based on high/low std, all trials')
xlabel(' delta F/F hit trials');
ylabel('delta F/F miss trials');


%%
AllAnimalsStimCluster = [];
AllAnimalsAllCells = [];
for r = 1:size(Rec,1)
    AllAnimalsStimCluster = [AllAnimalsStimCluster;mean(Activity(r).all(clusters(r).StimCluster,:,Performance(r).hit ),3) - mean(Activity(r).all(clusters(r).StimCluster,:,Performance(r).miss),3)]    ;
    AllAnimalsAllCells = [AllAnimalsAllCells;mean(Activity(r).all(:,:,Performance(r).hit ),3) - mean(Activity(r).all(:,:,Performance(r).miss),3)]    ;
end
Animal1DataM = mean(Activity(1).all(clusters(1).StimCluster,:,Performance(1).hit ),3) - mean(Activity(1).all(clusters(1).StimCluster,:,Perfolookrmance(1).miss),3);
Animal2DataM = mean(Activity(2).all(clusters(2).StimCluster,:,Performance(2).hit ),3) - mean(Activity(2).all(clusters(2).StimCluster,:,Performance(2).miss),3);
Animal3DataM = mean(Activity(3).all(clusters(3).StimCluster,:,Performance(3).hit ),3) - mean(Activity(3).all(clusters(3).StimCluster,:,Performance(3).miss),3);
Animal1Data = Activity(1).all(clusters(1).StimCluster,:,Stimuli{1,1}.TrialsPosition1 ) ;
Animal2Data = Activity(2).all(clusters(2).StimCluster,:,Stimuli{1,2}.TrialsPosition1) ;
Animal3Data = Activity(3).all(clusters(3).StimCluster,:,Stimuli{1,3}.TrialsPosition1) ;
Animal1LickData = Performance(1).LicksInFrame(:,Stimuli{1,1}.TrialsPosition1)
Animal2LickData = Performance(2).LicksInFrame(:,Stimuli{1,2}.TrialsPosition1)
Animal3LickData = Performance(3).LicksInFrame(:,Stimuli{1,3}.TrialsPosition1)
Animal1ReactionTime = []
Animal2ReactionTime = []
Animal3ReactionTime = []
for a = 1:size(Animal1LickData,2)
    Animal1ReactionTime = [Animal1ReactionTime, Animal1LickData(find(Animal1LickData(:,a) > 0,1),a)];
end
for a = 1:size(Animal2LickData,2)
    Animal2ReactionTime = [Animal2ReactionTime, Animal2LickData(find(Animal2LickData(:,a) > 0,1),a)];
end
for a = 1:size(Animal3LickData,2)
    Animal3ReactionTime = [Animal3ReactionTime, Animal3LickData(find(Animal3LickData(:,a) > 0,1),a)];
end
%%
AllAnimalsStimCluster = [];
AllAnimalsNonStimCluster = [];
AnimalReactionTimes = [];
timeReference = linspace(Time.Start,Time.Stop,Time.WindowLength);
ifi = 18.5/Time.WindowLength;
LickStartFrame = [];
PreLickModulationStimCluster = [];
PostLickModulationStimCluster = [];
PreLickModulationNonStimCluster = [];
PostLickModulationNonStimCluster = [];
PreRewardModulationStimCluster = [];
PostRewardModulationStimCluster = [];
PreRewardModulationNonStimCluster = [];
PostRewardModulationNonStimCluster = [];



for i = 1:size(Rec,1)
    AllAnimalsStimCluster = [AllAnimalsStimCluster; {Activity(i).all(clusters(i).StimCluster,:,Performance(i).hit) - mean(Activity(i).all(clusters(i).StimCluster,:,Performance(i).miss),3)}]    ;
    AllAnimalsNonStimCluster = [AllAnimalsNonStimCluster;{Activity(i).all(clusters(i).NonStimCluster,:,Performance(i).hit ) - mean(Activity(i).all(clusters(i).NonStimCluster,:,Performance(i).miss),3)}]    ;
    RecordingLickData = Performance(i).LicksInFrame(:,Performance(i).hit);
    FirstLick = [];
    StartFrame = [];
    PreL1 = [];
    PreL2 = [];
    PostL1 = [];
    PostL2 = [];
    PreR1 = [];
    PreR2 = [];
    PostR1 = [];
    PostR2 = [];

    for r = 1:size(RecordingLickData,2)
        temp = Performance(i).LicksInFrame(find(RecordingLickData(:,r) > 0 ,1),Performance(i).hit(r));
        if temp < 2
            FirstLick = [FirstLick;   temp];
            StartFrame = [StartFrame;find(FirstLick(end) < timeReference,1 )]
            PreL1 = cat(2,PreL1,AllAnimalsStimCluster{i,1}(:,StartFrame,r));
            PreL2 = cat(2,PreL2,AllAnimalsNonStimCluster{i,1}(:,StartFrame,r));
            PostL1 = cat(2,PostL1,AllAnimalsStimCluster{i,1}(:,StartFrame+2,r));
            PostL2 = cat(2,PostL2,AllAnimalsNonStimCluster{i,1}(:,StartFrame+2,r));
            PreR1 = cat(2,PreR1,AllAnimalsStimCluster{i,1}(:,StartFrame+3,r));
            PreR2 = cat(2,PreR2,AllAnimalsNonStimCluster{i,1}(:,StartFrame+3,r));
            PostR1 = cat(2,PostR1,AllAnimalsStimCluster{i,1}(:,StartFrame+6,r));
            PostR2 = cat(2,PostR2,AllAnimalsNonStimCluster{i,1}(:,StartFrame+5,r));
        end
    end
    AnimalReactionTimes = [AnimalReactionTimes;{FirstLick}];
    LickStartFrame= [LickStartFrame; {StartFrame}];
    PreLickModulationStimCluster = [PreLickModulationStimCluster;{PreL1}]
    PreLickModulationNonStimCluster = [PreLickModulationNonStimCluster;{PreL2}]
    PostLickModulationStimCluster = [PostLickModulationStimCluster;{PostL1}]
    PostLickModulationNonStimCluster = [PostLickModulationNonStimCluster;{PostL2}]
    PreRewardModulationStimCluster = [PreRewardModulationStimCluster;{PreR1}]
    PreRewardModulationNonStimCluster = [PreRewardModulationNonStimCluster;{PreR2}]
    PostRewardModulationStimCluster = [PostRewardModulationStimCluster;{PostR1}]
    PostRewardModulationNonStimCluster = [PostRewardModulationNonStimCluster;{PostR2}]
end

% figure;scatter(vertcat(PreLickModulationStimCluster{:}),vertcat(PostLickModulationStimCluster{:}))
% hold on;
% axis square
% hline1 = refline(1)
% hline1.LineWidth = 5
% hline1.Color ='r'
% xlabel(' Pre Lick Modulation');
% ylabel('Post Lick Modulation');

% figure;scatter(vertcat(PreLickModulationNonStimCluster{:}),vertcat(PostLickModulationNonStimCluster{:}))
% hold on;
% axis square
% hline1 = refline(1)
% hline1.LineWidth = 5
% hline1.Color ='r'
% xlabel(' Pre Lick Modulation');
% ylabel('Post Lick Modulation');

figure;%dscatter(vertcat(PreLickModulationStimCluster{:}),vertcat(PostLickModulationStimCluster{:}),'MSIZE',30)
hold on;
for a = 1:8
    scatter(mean(PreRewardModulationStimCluster{a, 1} ,2 ),mean(PostRewardModulationStimCluster{a, 1},2) ,30, [0 0 0],'filled')
end
axis square
hline1 = refline(1)
hline1.LineWidth = 5
hline1.Color ='r'
xlabel(' Pre Reward Modulation');
ylabel('Post Reward Modulation');
fontsize(24,'points')


%%
timeReference = linspace(Time.Start,Time.Stop,Time.WindowLength);
ifi = 18.5/Time.WindowLength;
PreExtraneousLick = [];
PostExtraneousLick = [];
for r = 1: size(Rec,1)
    for t = 1:size(Performance(r).LicksInFrame,2)
        extraneousLicks  = Performance(r).LicksInFrame(Performance(r).LicksInFrame(:,t) < -0.2 & Performance(r).LicksInFrame(:,t) > -7.9 ,t);
        if ~isempty(extraneousLicks)
            for e = 1:length(extraneousLicks)
                StartFrame = find(extraneousLicks(e) < timeReference,1 );
                PreExtraneousLick = [PreExtraneousLick;mean(Activity(r).all(:,StartFrame-2,t),1)];
                PostExtraneousLick = [PostExtraneousLick;mean(Activity(r).all(:,StartFrame+2,t),1)];
            end
        end
    end
end

figure;
hold on;
for a = 1:2
    scatter(PreExtraneousLick,PostExtraneousLick,30, [0 0 0],'filled')
end
axis square
hline1 = refline(1)
hline1.LineWidth = 5
hline1.Color ='r'
title(['PreLick / PostLick = ' num2str(mean(PreExtraneousLick./PostExtraneousLick))] )
xlabel(' dFF Pre Lick');
ylabel('dFF Post Lick');
fontsize(24,'points')
data = {};
data{1,1}(:,1) = PreExtraneousLick;
data{1,1}(:,2) = PostExtraneousLick;
names = ["Before Lick" , "After Lick"]
c =[params.colour(13,:)];
group_inx = repmat(1,size(data{1},1),1) ;
%
fig=figure('Position',[ 134         160        1311         675]);
% half-violin plots with white boxplots, jittered scatter and linkline
h = daviolinplot(data,'groups',group_inx,'colors',c,'box',3,...
    'boxcolor','w','scatter',2,'jitter',1,'scattercolor','same',...
    'scattersize',20,'scatteralpha',0.7,'linkline',1,'withinlines',1,...
    'xtlabels', names);
ylabel('dF/F');
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]); % make more space for the legend
set(h.sc,'MarkerEdgeColor','none');      % remove marker edge color
set(gca,'FontSize',24);
%%
data = {};
% data{1,1}(:,1) = reshape(PreLickModulationStimCluster{1},1,[]);
% data{1,1}(:,2) = reshape(PostLickModulationStimCluster{1},1,[]);
data{1,1}(:,1) = mean(PreLickModulationStimCluster{1},2);
data{1,1}(:,2) = mean(PostLickModulationStimCluster{1},2);
data{1,1}(:,3) = mean(PostRewardModulationStimCluster{1},2);

data{1,2}(:,1) = mean(PreLickModulationNonStimCluster{1},2);
data{1,2}(:,2) = mean(PostLickModulationNonStimCluster{1},2);
data{1,2}(:,3) = mean(PostRewardModulationNonStimCluster{1},2);

names = ["Before Lick" , "After Lick", "After Reward"]
c =[params.colour(13,:);params.colour(11,:)];
group_inx = [repmat(1,size(data{1},1),1);repmat(2,size(data{2},1),1)] ;
%%


fig=figure('Position',[ 134         160        1311         675]);
% half-violin plots with white boxplots, jittered scatter and linkline
h = daviolinplot(data(:,1),'groups',group_inx(:,1),'colors',c(1,:),'box',3,...
    'boxcolor','w','scatter',2,'jitter',1,'scattercolor','same',...
    'scattersize',20,'scatteralpha',0.7,'linkline',1,'withinlines',1,...
    'xtlabels', names);
ylabel('dF/F modulation');
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]); % make more space for the legend
set(h.sc,'MarkerEdgeColor','none');      % remove marker edge color
set(gca,'FontSize',24);


%%
fig=figure('Position',[ 134         160        1311         675]);
% half-violin plots with white boxplots, jittered scatter and linkline
h = daviolinplot(data,'groups',group_inx,'colors',c,'box',3,...
    'boxcolor','w','scatter',2,'jitter',1,'scattercolor','same',...
    'scattersize',20,'scatteralpha',0.7,'linkline',1,'withinlines',1,...
    'xtlabels', names);
ylabel('dF/F modulation');
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]); % make more space for the legend
set(h.sc,'MarkerEdgeColor','none');      % remove marker edge color
set(gca,'FontSize',24);


%%

figure;dscatter(vertcat(PreRewardModulationStimCluster{:}),vertcat(PostRewardModulationStimCluster{:}),'MSIZE',30)
hold on;
axis square
hline1 = refline(1)
hline1.LineWidth = 5
hline1.Color ='r'
% xlabel(' Pre Lick Modulation');
% ylabel('Post Lick Modulation');
fontsize(24,'points')
xlim([-15,25])
ylim([-15,25])


figure;dscatter(vertcat(PreRewardModulationNonStimCluster{:}),vertcat(PostRewardModulationNonStimCluster{:}),'MSIZE',30)
hold on;
axis square
hline1 = refline(1)
hline1.LineWidth = 5
hline1.Color ='r'
% xlabel(' Pre Lick Modulation');
% ylabel('Post Lick Modulation');
fontsize(24,'points')
xlim([-15,25])
ylim([-15,25])


%%
tl = tiledlayout(2,2,'TileSpacing', 'tight');

nexttile(1);
dscatter(vertcat(PreLickModulationStimCluster{:}),vertcat(PostLickModulationStimCluster{:}),'MSIZE',30)
hold on;
axis square
hline1 = refline(1)
hline1.LineWidth = 5
hline1.Color ='r'
xlabel(' Pre Lick Modulation');
ylabel('Post Lick Modulation');

nexttile(3);
dscatter(vertcat(PreLickModulationNonStimCluster{:}),vertcat(PostLickModulationNonStimCluster{:}),'MSIZE',30)
hold on;
axis square
hline1 = refline(1)
hline1.LineWidth = 5
hline1.Color ='r'
xlabel(' Pre Lick Modulation');
ylabel('Post Lick Modulation');

nexttile(2);
dscatter(vertcat(PreRewardModulationStimCluster{:}),vertcat(PostRewardModulationStimCluster{:}),'MSIZE',30)
hold on;
axis square
hline1 = refline(1)
hline1.LineWidth = 5
hline1.Color ='r'
xlabel(' Pre Reward Modulation');
ylabel('Post Reward Modulation');

nexttile(4);
dscatter(vertcat(PreRewardModulationNonStimCluster{:}),vertcat(PostRewardModulationNonStimCluster{:}),'MSIZE',30)
hold on;
axis square
hline1 = refline(1)
hline1.LineWidth = 5
hline1.Color ='r'
xlabel(' Pre Reward Modulation');
ylabel('Post Reward Modulation');


%%
figure;dscatter(vertcat(PreLickModulationStimCluster{:}),vertcat(PostRewardModulationStimCluster{:}),'MSIZE',30)
hold on;
axis square
hline1 = refline(1)
hline1.LineWidth = 5
hline1.Color ='r'
%%

[R1,ReactionTimeSortedIndex1] = sort(AnimalReactionTimes{1, 1});
[R2,ReactionTimeSortedIndex2] = sort(AnimalReactionTimes{2, 1});
[R3,ReactionTimeSortedIndex3] = sort(AnimalReactionTimes{3, 1});
% figure;t1 =tiledlayout(2,1);title(t1,'Recording 1');nexttile; imagesc(linspace(Time.Start,Time.Stop,Time.WindowLength),1:size(Animal1Data,3),squeeze(mean(Animal1Data(:,:,ReactionTimeSortedIndex1),1))');
% title('Average response of stimulus cluster (sorted by reaction Time)');nexttile;scatter(flip(R1),1:size(R1,2)); xlim([-8 10]);ylim([1 length(R1)]); title('First lick recorded');h = gca;h.YAxis.Visible = 'off';
% colormap(params.cmap);
% figure;t2 = tiledlayout(2,1);title(t2,'Recording 2');nexttile; imagesc(linspace(Time.Start,Time.Stop,Time.WindowLength),1:size(Animal2Data,3),squeeze(mean(Animal2Data(:,:,ReactionTimeSortedIndex2),1))');
% title('Average response of stimulus cluster (sorted by reaction Time)');nexttile;scatter(flip(R2),1:size(R2,2)); xlim([-8 10]);ylim([1 length(R2)]); title('First lick recorded');h = gca;h.YAxis.Visible = 'off';
% colormap(params.cmap);
% figure;t3 = tiledlayout(2,1);title(t3,'Recording 3');nexttile; imagesc(linspace(Time.Start,Time.Stop,Time.WindowLength),1:size(Animal3Data,3),squeeze(mean(Animal3Data(:,:,ReactionTimeSortedIndex3),1))');
% title('Average response of stimulus cluster (sorted by reaction Time)');nexttile;scatter(flip(R3),1:size(R3,2)); xlim([-8 10]); ylim([1 length(R3)]);title('First lick recorded');h = gca;h.YAxis.Visible = 'off';
% colormap(params.cmap);

figure;title('Recording 1');nexttile; imagesc(linspace(Time.Start,Time.Stop,Time.WindowLength),1:size(Animal1Data,3),squeeze(mean(Animal1Data(:,:,ReactionTimeSortedIndex1),1))');
title('Average response of stimulus cluster (sorted by reaction Time)');hold on;scatter(flip(R1), size(R1,2):-1:1,'r'); xlim([-8 10]);ylim([1 length(R1)]);h = gca;h.YAxis.Visible = 'off';
colormap(params.cmap);
figure;title('Recording 2');nexttile; imagesc(linspace(Time.Start,Time.Stop,Time.WindowLength),1:size(Animal2Data,3),squeeze(mean(Animal2Data(:,:,ReactionTimeSortedIndex2),1))');
title('Average response of stimulus cluster (sorted by reaction Time)');hold on;scatter(flip(R2), size(R2,2):-1:1,'r'); xlim([-8 10]);ylim([1 length(R2)]);h = gca;h.YAxis.Visible = 'off';
colormap(params.cmap);
figure;title('Recording 3');nexttile; imagesc(linspace(Time.Start,Time.Stop,Time.WindowLength),1:size(Animal3Data,3),squeeze(mean(Animal3Data(:,:,ReactionTimeSortedIndex3),1))');
title('Average response of stimulus cluster (sorted by reaction Time)');hold on;scatter(flip(R3), size(R3,2):-1:1,'r'); xlim([-8 10]); ylim([1 length(R3)]);h = gca;h.YAxis.Visible = 'off';
colormap(params.cmap);

%% Lick Triggered Average
timeReference = linspace(Time.Start,Time.Stop,Time.WindowLength);
ifi = 18.5/Time.WindowLength;
LTAstartFrame = [];
LTAWindow = [];
for r = 1:size(Rec,1)
    SF = [];
    for t = 1:size(AllAnimalsStimCluster{r, 1},3)
        SF = [SF; find(AnimalReactionTimes{r, 1}(t) < timeReference,1)];
    end
    LTAstartFrame = [LTAstartFrame;{SF}]
    LTAWindow = [LTAWindow; {Time.WindowLength - max(SF)}];
end
LTAWindowAll = min(ceil(vertcat(LTAWindow{:})/2));
LTAallCells = [];
LTA= [];
LTAincludedReactionTime = [];
for r = 1:size(Rec,1)
    LTAtemp = [];
    for t = 1:size(AllAnimalsStimCluster{r, 1},3)
        if AnimalReactionTimes{r, 1}(t) >1
            LTAtemps =  cat(3,LTAtemp, AllAnimalsStimCluster{r, 1}(:,LTAstartFrame{r,1}(t)-LTAWindowAll : LTAstartFrame{r,1}(t) + LTAWindowAll,t));
            LTA = [LTA; mean(AllAnimalsStimCluster{r, 1}(:,LTAstartFrame{r,1}(t)-LTAWindowAll : LTAstartFrame{r,1}(t) + LTAWindowAll,t),1)];
            LTAincludedReactionTime = [ LTAincludedReactionTime, AnimalReactionTimes{r, 1}(t) ];
        end
    end
    LTAallCells = [LTAallCells;{LTAtemp}]
end
[R1,R2] = sort(LTAincludedReactionTime)
figure;imagesc(linspace(-LTAWindowAll*ifi,LTAWindowAll*ifi, size(LTA,2)),1:size(LTA,1),LTA(R2,:));
colormap(params.cmap)
title('First lick triggered response (trials sorted by reaction time)');
xlabel('time (0 = first lick)');ylabel('Trial number');
xlim([-2 2])
xline(0,'-','First Lick');
xline(0.35,'-','Reward delivered');
fontsize(24,'points')

LTAM = mean(LTA,1);
figure;plot(linspace(-LTAWindowAll*ifi,LTAWindowAll*ifi, size(LTAM,2)),LTAM,'LineWidth',4);
title(' First Lick Triggered Average');
xlim([-2 2])
xline(0,'-','First Lick');
xline(0.35,'-','Reward delivered');
fontsize(24,'points')


%%
timeReference = linspace(Time.Start,Time.Stop,Time.WindowLength);
ifi = 18.5/Time.WindowLength;
LTAstartFrames = [];
LTAincludedTrials1 = find(Animal1ReactionTime > 0.5 & Animal1ReactionTime < 1);
LTAincludedTrials2 = find(Animal2ReactionTime > 0.5 & Animal2ReactionTime < 1);
LTAincludedTrials3 = find(Animal3ReactionTime > 0.5 & Animal3ReactionTime < 1);
for t = LTAincludedTrials1
    LTAstartFrames = [LTAstartFrames, find(Animal1ReactionTime(t) < timeReference,1 )];
end
for t = LTAincludedTrials2
    LTAstartFrames = [LTAstartFrames, find(Animal2ReactionTime(t) < timeReference,1 )];
end
for t = LTAincludedTrials3
    LTAstartFrames = [LTAstartFrames, find(Animal3ReactionTime(t) < timeReference,1 )];
end
LTAwindow = floor((Time.WindowLength - max(LTAstartFrames))/3);
LTAall = [];
LTAallTrials = [];

for t =  1:size(LTAincludedTrials1,2)
    % LTAall =  cat(3,LTAall, Animal1Data(:,LTAstartFrames(t)-LTAwindow : LTAstartFrames(t) + LTAwindow,LTAincludedTrials1(t)));
    LTAallTrials = [LTAallTrials; mean(Animal1Data(:,LTAstartFrames(t)-LTAwindow : LTAstartFrames(t) + LTAwindow,LTAincludedTrials1(t)),1)];
end
add = size(LTAincludedTrials1,2)
for t =  1:size(LTAincludedTrials2,2)
    % LTAall =  cat(3,LTAall, Animal2Data(:,LTAstartFrames(t+add)-LTAwindow : LTAstartFrames(t+add) + LTAwindow,LTAincludedTrials2(t)));
    LTAallTrials = [LTAallTrials; mean(Animal2Data(:,LTAstartFrames(t+add)-LTAwindow : LTAstartFrames(t+add) + LTAwindow,LTAincludedTrials2(t)),1)];
end
add = add+size(LTAincludedTrials2,2)
for t =  1:size(LTAincludedTrials3,2)
    % LTAall =  cat(3,LTAall, Animal3Data(:,LTAstartFrames(t+add)-LTAwindow : LTAstartFrames(t+add) + LTAwindow,LTAincludedTrials3(t)));
    LTAallTrials = [LTAallTrials; mean(Animal3Data(:,LTAstartFrames(t+add)-LTAwindow : LTAstartFrames(t+add) + LTAwindow,LTAincludedTrials3(t)),1)];
end

LTA_1to2 =  mean(LTAallTrials,1);
[ReactionTimesAllSorted,ReactionTimesAllSortedIndex] = sort([Animal1ReactionTime(LTAincludedTrials1),Animal2ReactionTime(LTAincludedTrials2),Animal3ReactionTime(LTAincludedTrials3)]);


fig=figure('Position',[ 134         160        1311         675]);
plot(linspace(-LTAwindow*ifi,LTAwindow*ifi, size(LTA_1to2,2)),LTA_1to2,'LineWidth',4,'Color','k');
xlim([-3 3])
title('')
xlabel('Time to first lick')
ylabel('dF/F')


fig=figure('Position',[ 134         160        1311         675]);
imagesc(linspace(-LTAwindow*ifi,LTAwindow*ifi, size(LTAallTrials,2)), 1:size(LTAallTrials,1), LTAallTrials(ReactionTimesAllSortedIndex,:));
colormap(params.cmap)

fig=figure('Position',[ 134         160        1311         675]);
title('Lick triggered averages')
f = tiledlayout(3,1);
ax1 = nexttile([2 1]);
imagesc(linspace(-LTAwindow*ifi,LTAwindow*ifi, size(LTA_1to2,2)),1: size(LTAallTrials,1), LTAallTrials(ReactionTimesAllSortedIndex,:));
colormap(params.cmap)
colorbar
ylabel('Trials')
xlabel('Time to first lick')
xlim([-3 3])
title('Lick triggered averages')
ax2 = nexttile;
plot(linspace(-LTAwindow*ifi,LTAwindow*ifi, size(LTA_1to2,2)),LTA_1to2,'LineWidth',4,'Color','k');
xlim([-3 3])
xlabel('Time to first lick')
ylabel('dF/F')
fontsize(24,'points')





%% Lick Trigerred Modulation
timeReference = linspace(Time.Start,Time.Stop,Time.WindowLength);
ifi = 18.5/Time.WindowLength;
LTMstartFrames = [];
LTMincludedTrials1 = find(Animal1ReactionTime > 0.5 & Animal1ReactionTime < 1);
LTMincludedTrials2 = find(Animal2ReactionTime > 0.5 & Animal2ReactionTime < 1);
LTMincludedTrials3 = find(Animal3ReactionTime > 0.5 & Animal3ReactionTime < 1);
for t = LTMincludedTrials1
    LTMstartFrames = [LTMstartFrames, find(Animal1ReactionTime(t) < timeReference,1 )];
end
for t = LTMincludedTrials2
    LTMstartFrames = [LTMstartFrames, find(Animal2ReactionTime(t) < timeReference,1 )];
end
for t = LTMincludedTrials3
    LTMstartFrames = [LTMstartFrames, find(Animal3ReactionTime(t) < timeReference,1 )];
end
LTMwindow = floor((Time.WindowLength - max(LTMstartFrames))/2);
LTMall = [];
LTMallTrials = [];
[a,MissTrial] = ismember(Performance(1).miss , Stimuli{1, 1}.TrialsPosition1 )
MissAverage = mean(Animal1Data(:,:,MissTrial),[1,3])
[a,MissTrial] = ismember(Performance(2).miss , Stimuli{1, 2}.TrialsPosition1 )
MissAverage =[MissAverage; mean(Animal2Data(:,:,MissTrial),[1,3])]
[a,MissTrial] = ismember(Performance(3).miss , Stimuli{1, 3}.TrialsPosition1 )
MissAverage =[MissAverage; mean(Animal3Data(:,:,MissTrial),[1,3])]

for t =  1:size(LTMincludedTrials1,2)
    LTMallTrials = [LTMallTrials; mean(Animal1Data(:,LTMstartFrames(t)-LTMwindow : LTMstartFrames(t) + LTMwindow,LTMincludedTrials1(t)),1) - MissAverage(1,LTMstartFrames(t)-LTMwindow : LTMstartFrames(t) + LTMwindow)];
end
add = size(LTMincludedTrials1,2)
for t =  1:size(LTMincludedTrials2,2)
    LTMallTrials = [LTMallTrials; mean(Animal2Data(:,LTMstartFrames(t+add)-LTMwindow : LTMstartFrames(t+add) + LTMwindow,LTMincludedTrials2(t)),1) - MissAverage(2,LTMstartFrames(t)-LTMwindow : LTMstartFrames(t) + LTMwindow)];
end
add = add+size(LTMincludedTrials2,2)
for t =  1:size(LTMincludedTrials3,2)
    LTMallTrials = [LTMallTrials; mean(Animal3Data(:,LTMstartFrames(t+add)-LTMwindow: LTMstartFrames(t+add) +LTMwindow,LTMincludedTrials3(t)),1) - MissAverage(3,LTMstartFrames(t)-LTMwindow: LTMstartFrames(t) +LTMwindow)];
end

LTM =  mean(LTMallTrials,1);
[ReactionTimesAllSorted,ReactionTimesAllSortedIndex] = sort([Animal1ReactionTime(LTMincludedTrials1),Animal2ReactionTime(LTMincludedTrials2),Animal3ReactionTime(LTMincludedTrials3)]);


% fig=figure('Position',[ 134         160        1311         675]);
% plot(linspace(-LTMwindow*ifi,LTMwindow*ifi, size(LTM,2)),LTM,'LineWidth',4,'Color','k');
% xlim([-1 1])
% title('')
% xlabel('Time to first lick')
% ylabel('dF/F')
%
%
% fig=figure('Position',[ 134         160        1311         675]);
% imagesc(linspace(-LTMwindow*ifi,LTMwindow*ifi, size(LTM,2)), 1:size(LTMallTrials,1), LTMallTrials(ReactionTimesAllSortedIndex,:));
% colormap(params.cmap)

fig=figure('Position',[ 134         160        1311         675]);
title('Lick triggered modulation')
f = tiledlayout(3,1);
ax1 = nexttile([2 1]);
imagesc(linspace(-LTMwindow*ifi,LTMwindow*ifi, size(LTM,2)), 1:size(LTMallTrials,1), LTMallTrials(ReactionTimesAllSortedIndex,:));
colormap(params.cmap)
colorbar
clim([-2 2])
ylabel('Trials')
xlabel('Time to first lick')
xlim([-4 4])
title('Lick triggered modulation')
ax2 = nexttile;
plot(linspace(-LTMwindow*ifi,LTMwindow*ifi, size(LTM,2)),LTM,'LineWidth',4,'Color','k');
xlim([-4 4])
xlabel('Time to first lick')
ylabel('dF/F')
fontsize(24,'points')





%%

figure;
aa = tiledlayout(1,2);
nexttile;
imagesc(linspace(-LTAwindow1*ifi,LTAwindow1*ifi, size(LTA1allTrials,2)),1:size(LTA1allTrials,1),LTA1allTrials(ReactionTimeSortedIndex1,:));
colormap(params.cmap)
hold on;scatter(-flip(R1),size(R1,2):-1:1,300,'r','filled')
title(aa,'Lick triggered response','fontsize',24);
xlabel('time (0 = first lick)');ylabel('Trial number');
nexttile;
imagesc(linspace(-LTAwindow3*ifi,LTAwindow3*ifi, size(LTA3allTrials,2)),1:size(LTA3allTrials,1),LTA3allTrials(ReactionTimeSortedIndex3,:));
colormap(params.cmap)
hold on;scatter(-flip(R3),size(R3,2):-1:1,300,'r','filled')
xlabel('time (0 = first lick)');ylabel('Trial number');
fontsize(24,'points')




%% influence of previous trial
PreviousHit =[]
PreviousMiss =[]
temp = find(Animal1ReactionTime < 2 == 1);
temp(temp == 1) = [];
Performance(1).ph = temp - 1
Performance(1).pm = setdiff(1:size(temp,2),temp)
temp = find(Animal2ReactionTime < 2 == 1);
temp(temp == 1) = [];
Performance(2).ph = temp - 1
Performance(2).pm = setdiff(1:size(temp,2),temp)
temp = find(Animal3ReactionTime < 2 == 1);
temp(temp == 1) = [];
Performance(3).ph = temp - 1
Performance(3).pm = setdiff(1:size(temp,2),temp)


for r = 1:2%size(Rec,1)
    PreviousHit = [PreviousHit;mean(Activity(r).all(clusters(r).StimCluster,:,Performance(r).ph),3)];
    PreviousMiss = [PreviousMiss;mean(Activity(r).all(clusters(r).StimCluster,:,Performance(r).pm),3)];
end

nrSamples = 2;
cMap = lines(nrSamples);
figure;
hold on;
aline(2) = stdshade(PreviousMiss,0.15,params.colour(3,:),linspace(Time.Start,Time.Stop,Time.WindowLength));
aline(1) = stdshade(PreviousHit,0.15,params.colour(1,:),linspace(Time.Start,Time.Stop,Time.WindowLength));
aline(2).LineWidth  = 3
aline(1).LineWidth  = 3
legend(aline,'previous Hit','previous Miss');
axis square
xlabel('Time in seconds');
ylabel('delta F/F');
title(params.titleName);
subtitle(join([params.AnimalName ,',' ,params.TrainingState]))
xlim([-2 4])
fontsize(24,"points")

%%
saveThis = 0;


for i = 1:length(Performance(1).hit)
    aa = figure;imagesc(linspace(Time.Start,Time.Stop,Time.WindowLength),1:size(Animal1Data,1),Activity(1).all(clusters(1).StimCluster,:,Performance(1).hit(i)) - mean(Activity(1).all(clusters(1).StimCluster,:,Performance(1).miss),3));
    hold on;line([Animal1ReactionTime(find(Performance(1).hit(i) == Stimuli{1,1}.TrialsPosition1)) Animal1ReactionTime(find(Performance(1).hit(i) == Stimuli{1,1}.TrialsPosition1))],[1 size(LTA1all,1)],'LineWidth',2,'color','k')
    xlabel('time'); ylabel('neuron ID');
    title('Single hit trial - average miss trial response', 'fontsize',24)
    % clim([-1 2])
    colormap(params.cmap)
    if saveThis
        if i == 1
            exportgraphics(aa,[paths.outputdir '\SingleTrialMinusMiss.pdf'])
        else
            exportgraphics(aa,[paths.outputdir '\SingleTrialMinusMiss.pdf'],'Append',true)
        end
    end
end

%%
figure;
aa = tiledlayout(2,1);

nexttile;
imagesc(linspace(-2,4,61),1:size(Animal1Data,1),Activity(1).all(clusters(1).StimCluster,62:121,Performance(1).hit(9)) - mean(Activity(1).all(clusters(1).StimCluster,62:121,Performance(1).miss),3));
hold on;line([Animal1ReactionTime(find(Performance(1).hit(9) == Stimuli{1,1}.TrialsPosition1)) Animal1ReactionTime(find(Performance(1).hit(9) == Stimuli{1,1}.TrialsPosition1))],[1 size(LTA1all,1)],'LineWidth',2,'color','k')
xlabel('time'); ylabel('neuron ID');
nexttile;
imagesc(linspace(-2,4,61),1:size(Animal1Data,1),Activity(1).all(clusters(1).StimCluster,62:121,Performance(1).hit(7)) - mean(Activity(1).all(clusters(1).StimCluster,62:121,Performance(1).miss),3));
hold on;line([Animal1ReactionTime(find(Performance(1).hit(7) == Stimuli{1,1}.TrialsPosition1)) Animal1ReactionTime(find(Performance(1).hit(7) == Stimuli{1,1}.TrialsPosition1))],[1 size(LTA1all,1)],'LineWidth',2,'color','k')
xlabel('time'); ylabel('neuron ID');
% clim([-1 2])
colormap(params.cmap)
title(aa,'Single hit trial - average miss trial response', 'fontsize',24)

%%
figure;
imagesc(linspace(-2,4,61),1:size(Animal1Data,1),Activity(1).all(clusters(1).StimCluster,62:121,Performance(1).hit(7)) - mean(Activity(1).all(clusters(1).StimCluster,62:121,Performance(1).miss),3));
hold on;line([Animal1ReactionTime(find(Performance(1).hit(7) == Stimuli{1,1}.TrialsPosition1)) Animal1ReactionTime(find(Performance(1).hit(7) == Stimuli{1,1}.TrialsPosition1))],[1 size(LTA1all,1)],'LineWidth',2,'color','k')
xlabel('Time to stimulus');
ylabel('neuron ID');
% clim([-1 2])
colormap(params.cmap)
clim([-10 10])
title('Single hit trial - average miss trial response', 'fontsize',24)
fontsize(24,'points')



