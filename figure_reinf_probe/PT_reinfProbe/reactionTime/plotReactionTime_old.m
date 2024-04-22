function plotReactionTime(mouseMega)
biasThreshold = 0.4;
modelBias = mouseMega.getProp('behav','field','modelBias');
stateFlag = cellfun(@(x)(abs(x)>biasThreshold),modelBias,'UniformOutput',false);
reactionTime = mouseMega.getProp('behav','field','reactionTime');
reactionTime = cellfun(@(x)(fillNan(x,x>2.5)),reactionTime,'UniformOutput',false);
reactionTime = cellfun(@(x)(fillNan(x,x<0)),reactionTime,'UniformOutput',false);

%-------------------- BY Day--------------------------------------
day = mouseMega.getProp('behav','field','day');
trialNumThre = 30;
unbiasRT_day = cell(1,length(day)); biasRT_day = cell(1,length(day));
for i = 1:length(day)
    tempRT = smoothdata(reactionTime{i},'movmean',50); 
    figure; plot(tempRT,'LineWidth',2); hold on;
    tempRT_bias = tempRT; tempRT_bias(~stateFlag{i}) = nan;
    plot(tempRT_bias,'LineWidth',2); hold on;
    dayLim = max(day{i});
    %dayLim = 10;
    for j = 1:dayLim
        allTrialThisDay = sum(day{i}==j);
        biasTrialThisDay = sum(stateFlag{i} & (day{i}==j));
        if biasTrialThisDay>trialNumThre && (allTrialThisDay-biasTrialThisDay)>trialNumThre
            tempRT = reactionTime{i}(stateFlag{i} & (day{i}==j));
            biasRT_day{i} = cat(1,biasRT_day{i},nanmean(tempRT));
            tempRT = reactionTime{i}(~stateFlag{i} & (day{i}==j));
            unbiasRT_day{i} = cat(1,unbiasRT_day{i},nanmean(tempRT));
        end
    end
end
figure;  plotRT_bySession(biasRT_day,unbiasRT_day)
figure;  plotRT_byAni(biasRT_day,unbiasRT_day)

%------------- all trials
trialLim = 1400;
stateFlag = cellfun(@(x)(x(1:trialLim)),stateFlag,'UniformOutput',false);
reactionTime = cellfun(@(x)(x(1:trialLim)),reactionTime,'UniformOutput',false);


biasRT = cellfun(@(x,y)(x(y~=0)),reactionTime, stateFlag,'UniformOutput',false);
biasRT = cellfun(@nanmean,biasRT,'UniformOutput',true);

unbiasRT = cellfun(@(x,y)(x(y==0)),reactionTime, stateFlag,'UniformOutput',false);
unbiasRT = cellfun(@nanmean,unbiasRT,'UniformOutput',true);
    
figure; 
plot([0 1.5],[0 1.5],'Color',[0.8 0.8 0.8],'LineWidth',2);hold on; 
scatter(biasRT,unbiasRT,30,matlabColors(1),'filled'); 




end

function mat = fillNan(mat,flag)
    mat(flag) = nan;    
end

function plotRT_bySession(biasRT_day,unbiasRT_day)
    biasRT_day_allSession = fn_cell2mat(biasRT_day,1);
    unbiasRT_day_allSession = fn_cell2mat(unbiasRT_day,1);
    plot([0 1.5],[0 1.5],'Color',[0.8 0.8 0.8],'LineWidth',2);hold on; 
    scatter(biasRT_day_allSession,unbiasRT_day_allSession,30,matlabColors(1),'filled'); 
    xlabel('biased RT'); ylabel('unbiased RT')
end

function plotRT_byAni(biasRT_day,unbiasRT_day)
    biasRT_day = cellfun(@mean,biasRT_day);
    unbiasRT_day = cellfun(@mean,unbiasRT_day);
    plot([0 1.5],[0 1.5],'Color',[0.8 0.8 0.8],'LineWidth',2);hold on; 
    scatter(biasRT_day,unbiasRT_day,30,matlabColors(1),'filled'); 
    xlabel('biased RT'); ylabel('unbiased RT')
end