function LocomotionComparisionAndCorrelation(ActivityData, Behav,WindowLength,Rec, Performance, LocomotionCal, params)


Naive = []; Beginner = []; Expert = []; NoSpout = [];

for i = 1:length(Behav.Naive)
 Naive = [Naive,mean(Behav.Naive(i).allPupilSize,1)];
end
for i = 1:length(Behav.Beginner)
 Beginner = [Beginner,mean(Behav.Beginner(i).allPupilSize,1)];
end
for i = 1:length(Behav.Expert)
Expert = [Expert,mean(Behav.Expert(i).allPupilSize,1)]  ;
end
for i = [1,2,4]%length(Behav.NoSpout)
NoSpout = [NoSpout, mean(Behav.NoSpout(i).allPupilSize,1)];
end
data = [Naive,Beginner,Expert,NoSpout];
data = data / 1000;

names = ["Naive", "Beginner", "Expert","No Spout"];
c =[params.Naive;params.Beginner; params.Expert; params.NoSpout];

group_inx = [repmat(1,length(Naive),1);repmat(12,length(Beginner),1);repmat(3,length(Expert),1);repmat(4,length(NoSpout),1) ]';
figure
% half-violin plots with white boxplots, jittered scatter and linkline 
h = daviolinplot(data,'groups',group_inx,'colors',c,'box',2,'boxcolors' ,'same','xtlabels', names,'violin','full'); 
ylabel('Pupil Area [mm^2]');
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]); % make more space for the legend
set(gca,'FontSize',24);   
%%
Naive = []; Beginner = []; Expert = []; NoSpout = [];
Naive_mean = []; B_mean = []; E_mean = []; NoSpout_mean = [];

    for  a = 1:size(LocomotionCal.Naive,2)
   TrialLength = floor(length(LocomotionCal.Naive(a).Forward_mm) / 150);
        TrialMovement =   1:TrialLength:length(LocomotionCal.Naive(a).Forward_mm);
                Naive = [Naive;LocomotionCal.Naive(a).Forward_mm];        
         for t = 1:150
                      Naive_mean = [Naive_mean,mean(LocomotionCal.Naive(a).Forward_mm(TrialMovement(t):TrialMovement(t+1)) ,1)];
         end
    end

for i = 1:2%length(Behav.Beginner)
    % if not(isempty(LocomotionCal.Beginner{1, i}(1).allForward_mm ))
    for  a = 1:size(LocomotionCal.Beginner{1, i},2)
        B_mean = [B_mean,mean(LocomotionCal.Beginner{1, i}(a).Forward_mm ,1)];
        Beginner = [Beginner;LocomotionCal.Beginner{1, i}(a).Forward_mm];
    end
end

for i = 1:3%length(Behav.Expert)
    for  a = 1:size(LocomotionCal.Expert{1, i},2)
         E_mean = [E_mean,mean(LocomotionCal.Expert{1, i}(a).Forward_mm ,1)]  ;
                Expert = [Expert;LocomotionCal.Expert{1, i}(a).Forward_mm ]  ;

    end
end

  for  a = 1:size(LocomotionCal.NoSpout,2)
   TrialLength = floor(length(LocomotionCal.NoSpout(a).Forward_mm) / 150);
        TrialMovement =   1:TrialLength:length(LocomotionCal.NoSpout(a).Forward_mm);
                NoSpout = [NoSpout;LocomotionCal.NoSpout(a).Forward_mm];        
         for t = 1:150
                      NoSpout_mean = [NoSpout_mean,mean(LocomotionCal.NoSpout(a).Forward_mm(TrialMovement(t):TrialMovement(t+1)) ,1)];
         end
    end

data = [Naive; Beginner;Expert;NoSpout];
data = data * TrialLength / 12;
 data(find(data < 0)) = NaN;
  data(find(data > 25)) = 1;
data_mean = [Naive_mean, B_mean,E_mean,NoSpout_mean];
data_mean = data_mean * TrialLength / 12;
 data_mean(find(data_mean < 0)) = NaN;
  data_mean(find(data_mean > 25)) = 1;

names = ["Naive", "Beginner", "Expert","No Spout"];
c =[params.Naive;params.Beginner; params.Expert; params.NoSpout];
group_inx = [repmat(1,length(Naive),1);repmat(2,length(Beginner),1);repmat(3,length(Expert),1);repmat(4,length(NoSpout),1) ]';
    fig=figure('Position',[ 134         160        1311         675]);
% half-violin plots with white boxplots, jittered scatter and linkline 
h = daviolinplot(data,'groups',group_inx,'colors',c,'box',2,'boxcolors' ,'same','xtlabels', names,'violin','full'); 
ylabel('Speed [mm/sec]');
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]); % make more space for the legend
set(gca,'FontSize',24);   
ylim([-1 15])

% group_inx_mean = [repmat(1,length(Naive_mean),1);repmat(2,length(B_mean),1);repmat(3,length(E_mean),1);repmat(4,length(NoSpout_mean),1) ]';
%     fig=figure('Position',[ 134         160        1311         675]);
% % half-violin plots with white boxplots, jittered scatter and linkline 
% h = daviolinplot(data_mean,'groups',group_inx,'colors',c,'box',2,'boxcolors' ,'same','xtlabels', names,'violin','full'); 
% ylabel('Speed [mm/sec]');
% xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]); % make more space for the legend
% set(gca,'FontSize',24);   
% ylim([-1 15])
%%
figure
 violins = violinplot(group_inx,data)
ylim([-1 15])
xticklabels({"Naive", "Beginner", "Expert","No Spout"})
 % 
 % figure
 %  violins = violinplot(group_inx_mean, data_mean)
 % ylim([-1 15])
 % xticklabels({"Naive", "Beginner", "Expert","No Spout"})
%%
end