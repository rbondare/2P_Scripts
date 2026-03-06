%% Lick correlations
function lickCorrelations(ActivityData, Entropy, WindowLength, params, Rec, Performance)
% LICKCORRELATIONS  Analyse pre- and post-lick activity/entropy and produce
% correlation and distribution plots.
%   • Figure 1 : Pre vs Post dF/F scatter   – Wilcoxon signed‑rank (paired)
%   • Figure 2 : Pre vs Post Entropy scatter – Wilcoxon signed‑rank (paired)
%   • Figure 3 : Histograms of Post/Pre ratios (separate tiles) – sign‑rank
%                against a median of 1
%   • Figure 4 : Pre value vs Change scatter (dF/F & Entropy) – Spearman ρ
%                with p‑value; also signed‑rank of change ≠ 0

%% ------------------------- Data Accumulation ------------------------ %%
% Time axis (–8 s to +10.5 s around lick)
timeReference   = linspace(-8,10.5,WindowLength.All);
ifi             = 18.5/WindowLength.All; %#ok<NASGU>

PreExtraneousLickActivity  = [];
PostExtraneousLickActivity = [];
PreExtraneousLickEntropy   = [];
PostExtraneousLickEntropy  = [];

for r = 1:size(Rec.Expert,1)
    for t = 1:size(Performance.Expert(r).LicksInFrame,2)
        extraneousLicks = Performance.Expert(r).LicksInFrame( ...
            Performance.Expert(r).LicksInFrame(:,t) < -0.2 & ...
            Performance.Expert(r).LicksInFrame(:,t) > -7.9 , t);

        for e = 1:numel(extraneousLicks)
            StartFrame = find(extraneousLicks(e) < timeReference, 1);
            if StartFrame > 2 && StartFrame < (size(ActivityData.Expert(r).all,2)-2)
                PreExtraneousLickActivity  = [PreExtraneousLickActivity;  mean(ActivityData.Expert(r).all(:,StartFrame-2,t),1)]; %#ok<AGROW>
                PostExtraneousLickActivity = [PostExtraneousLickActivity; mean(ActivityData.Expert(r).all(:,StartFrame+2,t),1)]; %#ok<AGROW>
                PreExtraneousLickEntropy   = [PreExtraneousLickEntropy;   Entropy.Expert.AllNeurons.RecordingEntropy_Raw{r,1}(t,StartFrame-2)]; %#ok<AGROW>
                PostExtraneousLickEntropy  = [PostExtraneousLickEntropy;  Entropy.Expert.AllNeurons.RecordingEntropy_Raw{r,1}(t,StartFrame+2)]; %#ok<AGROW>
            end
        end
    end
end

%% Convenience variables
ActivityRatio  = PostExtraneousLickActivity ./ PreExtraneousLickActivity;
EntropyRatio   = PostExtraneousLickEntropy  ./ PreExtraneousLickEntropy;
ActivityChange = PostExtraneousLickActivity - PreExtraneousLickActivity;
EntropyChange  = PostExtraneousLickEntropy  - PreExtraneousLickEntropy;

%% ------------------------- Figure 1  dF/F scatter -------------------- %%
figure; hold on;
scatter(PreExtraneousLickActivity, PostExtraneousLickActivity, 30, [0 0 0], 'filled');
axis square;
ref = refline(1); ref.LineWidth = 5; ref.Color = 'r';

[pAct,~,statsAct] = signrank(PreExtraneousLickActivity, PostExtraneousLickActivity);
ratioAct = mean(PreExtraneousLickActivity ./ PostExtraneousLickActivity, 'omitnan');

ax = gca;
ax.FontSize = 20;

title({'Pre vs Post dF/F',sprintf('ratio = %.3f,  p = %.3g', ratioAct, pAct)});

%% Half‑violin (dF/F)
figure;
axesFont = 24;
daviolinplot({[PreExtraneousLickActivity, PostExtraneousLickActivity]}, ...
             'groups', ones(numel(PreExtraneousLickActivity),1), ...
             'colors', params.colour(13,:), 'box',3,'boxcolor','w', ...
             'scatter',2,'jitter',1,'scattercolor','same','scattersize',20, ...
             'scatteralpha',0.7,'linkline',0,'withinlines',0, ...
             'xtlabels',{"Before Lick","After Lick"});
set(gca,'FontSize',axesFont);
title(gca,sprintf('p = %.3g',pAct));

%% ------------------------- Figure 2  Entropy scatter ----------------- %%
figure; hold on;
scatter(PreExtraneousLickEntropy, PostExtraneousLickEntropy, 30, [0 0 0], 'filled');
axis square;
ref = refline(1); ref.LineWidth = 5; ref.Color = 'r';

[pEnt,~,statsEnt] = signrank(PreExtraneousLickEntropy, PostExtraneousLickEntropy);
ratioEnt = mean(PreExtraneousLickEntropy ./ PostExtraneousLickEntropy, 'omitnan');

ax = gca;
ax.FontSize = 24;

title(ax,{'Pre vs Post Entropy', sprintf('ratio = %.3f,  p = %.3g', ratioEnt, pEnt)});

%% Half‑violin (Entropy)
figure;
daviolinplot({[PreExtraneousLickEntropy, PostExtraneousLickEntropy]}, ...
             'groups', ones(numel(PreExtraneousLickEntropy),1), ...
             'colors', params.colour(13,:), 'box',3,'boxcolor','w', ...
             'scatter',2,'jitter',1,'scattercolor','same','scattersize',20, ...
             'scatteralpha',0.7,'linkline',1,'withinlines',0, ...
             'xtlabels',{"Before Lick","After Lick"});
set(gca,'FontSize',axesFont);
title(gca,sprintf('p = %.3g',pEnt));

%% ------------------------- Figure 3  Ratio Histograms --------------- %%
figure;
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% dF/F ratio tile
ax1 = nexttile;
histogram(ax1, ActivityRatio, 'Normalization','probability', ...
          'BinMethod','sturges', 'FaceColor', params.colour(7,:), ...
          'EdgeColor','none', 'FaceAlpha',0.7);
[pRatioAct] = signrank(ActivityRatio - 1);
ax1.FontSize = axesFont;
ax1.Box = 'on';
ax1.XLabel.String = 'Post / Pre dF/F';
ax1.YLabel.String = 'Probability';

title(ax1,{'dF/F Ratio',sprintf('p = %.3g', pRatioAct)});

% Entropy ratio tile
ax2 = nexttile;
histogram(ax2, EntropyRatio, 'Normalization','probability', ...
          'BinMethod','sturges', 'FaceColor', params.colour(9,:), ...
          'EdgeColor','none', 'FaceAlpha',0.7);
[pRatioEnt] = signrank(EntropyRatio - 1);
ax2.FontSize = axesFont;
ax2.Box = 'on';
ax2.XLabel.String = 'Post / Pre Entropy';
ax2.YLabel.String = 'Probability';

title(ax2,{'Entropy Ratio',sprintf('p = %.3g', pRatioEnt)});

%% ------------------------- Figure 4  Change vs Pre ------------------ %%
figure;

% dF/F change
ax3 = subplot(1,2,1);
scatter(ax3, PreExtraneousLickActivity, ActivityChange, 30, params.colour(7,:), 'filled');
lsline(ax3);
axis(ax3,'square'); grid(ax3,'on');
[rhoAct, pRhoAct] = corr(PreExtraneousLickActivity, ActivityChange, ...
                         'Type','Spearman','Rows','complete');
ax3.FontSize = axesFont;
ax3.XLabel.String = 'Pre dF/F';
ax3.YLabel.String = '\Delta dF/F (Post - Pre)';

title(ax3,{'Activity Change vs Pre',sprintf('rho = %.2f,  p = %.3g', rhoAct, pRhoAct)});

% Entropy change
ax4 = subplot(1,2,2);
scatter(ax4, PreExtraneousLickEntropy, EntropyChange, 30, params.colour(9,:), 'filled');
lsline(ax4);
axis(ax4,'square'); grid(ax4,'on');
[rhoEnt, pRhoEnt] = corr(PreExtraneousLickEntropy, EntropyChange, ...
                         'Type','Spearman','Rows','complete');
ax4.FontSize = axesFont;
ax4.XLabel.String = 'Pre Entropy';
ax4.YLabel.String = '\Delta Entropy (Post - Pre)';

title(ax4,{'Entropy Change vs Pre',sprintf('rho = %.2f,  p = %.3g', rhoEnt, pRhoEnt)});

end
