clear;
sep = '\';
dataPath = 'C:\Users\zzhu34\Documents\gitRep\octoBehavior\trialData';
mouse = 'zz068';
load([dataPath sep mouse '.mat']);

%% Get accuracy and action rate
trialSmoothBin = 30;
trialPerDay = [];
for i = 1:length(trialData.responseType)
    if ~isempty(trialData.responseType{i}); trialPerDay = [trialPerDay length(trialData.responseType{i})]; end 
end
cumTrialPerDay = cumsum(trialPerDay);

accuracy = cellfun(@(x)(double(x==1)),trialData.responseType,'UniformOutput',false);

accuracySmooth = cellfun(@(x)(smoothdata(x,'movmean',trialSmoothBin)),accuracy,'UniformOutput',false);
accuracySmooth = fn_cell2mat(accuracySmooth,1);

actionRate = cellfun(@(x)(double(x~=0)),trialData.responseType,'UniformOutput',false);
actionRateSmooth = cellfun(@(x)(smoothdata(x,'movmean',trialSmoothBin)),actionRate,'UniformOutput',false);
actionRateSmooth = fn_cell2mat(actionRateSmooth,1);

accuracySmooth(actionRateSmooth<0.4) = nan;
figure; hold on;
plot(accuracySmooth); hold on ; plot(accuracySmooth./actionRateSmooth);
for i = 1:length(trialPerDay)-1
    plot([cumTrialPerDay(i)+1 cumTrialPerDay(i)+1], [0 1],'Color',[0.8 0.8 0.8]);
end 
xlim([1 length(accuracySmooth)])
legend('acc.','acc. miss corr','Location','Best')

%day = 13;
%mdl = getPrevTrialRegression(trialData.stimulus{day},trialData.action{day},1,200,20);



function rSq = getPrevTrialRegression(stimulus,action,nTrialBack,binSize,windowSize)

stimulus(stimulus==2) = -1;
action(action==2) = -1;

if nargin<=3 % regression using all trials in the session
    rSq = getLogisticRegression(stimulus,action,nTrialBack);
    %{
    regressor = [];
    regressor = cat(2,regressor,stimulus); % regressor 1: current trial type

    reward = stimulus .* action; 
    posReward = double(reward==1); negReward = double(reward==-1);
    action_posReward = action .* posReward; action_negReward = action .* negReward;

    for i = 1:nTrialBack
        regressor = cat(2,regressor, circshift(stimulus,i)); % regressor 2: previous trial stimulus
        regressor = cat(2,regressor, circshift(action,i)); % regressor 3: previous trial actions
        %regressor = cat(2,regressor, circshift(action_posReward,i)); % regressor 4: previous pos reward X action association
        %regressor = cat(2,regressor, circshift(action_negReward,i)); % regressor 5: previous neg reward X action association
        regressor = cat(2,regressor, stimulus .* circshift(action_posReward,i)); % regressor 4: previous pos reward X action X current stimulus
        regressor = cat(2,regressor, stimulus .* circshift(action_negReward,i)); % regressor 5: previous pos reward X action X current stimulus
    end

    targVar = action(nTrialBack+1:end); regressor = regressor(nTrialBack+1:end,:);
    missFlag = (targVar ==0); targVar(missFlag) = []; regressor(missFlag,:) = [];
    regressor = zscore(regressor,1);
    targVar(targVar==1) = 0; targVar(targVar==-1) = 1;
    
    mdl = fitglm(regressor,targVar,'Distribution','binomial','Link','logit');
    %}

else % regression using blocks of trials, with a certain sliding window
    blockStart = 1:windowSize:(length(stimulus)-binSize); blockEnd = blockStart + binSize;
    for i = 1:length(blockStart)
        tempStim = stimulus(blockStart(i):blockEnd(i));
        tempAction = action(blockStart(i):blockEnd(i));
        mdl = getLogisticRegression(tempStim,tempAction,nTrialBack);
        rSq = mdl;
    end
    
end




end

function mdl = getLogisticRegression(stimulus,action,nTrialBack)
regressor = [];
regressor = cat(2,regressor,stimulus); % regressor 1: current trial type

reward = stimulus .* action; 
posReward = double(reward==1); negReward = double(reward==-1);
action_posReward = action .* posReward; action_negReward = action .* negReward;

for i = 1:nTrialBack
    regressor = cat(2,regressor, circshift(stimulus,i)); % regressor 2: previous trial stimulus
    regressor = cat(2,regressor, circshift(action,i)); % regressor 3: previous trial actions
    %regressor = cat(2,regressor, circshift(action_posReward,i)); % regressor 4: previous pos reward X action association
    %regressor = cat(2,regressor, circshift(action_negReward,i)); % regressor 5: previous neg reward X action association
    regressor = cat(2,regressor, stimulus .* circshift(action_posReward,i)); % regressor 4: previous pos reward X action X current stimulus
    regressor = cat(2,regressor, stimulus .* circshift(action_negReward,i)); % regressor 5: previous pos reward X action X current stimulus
end

targVar = action(nTrialBack+1:end); regressor = regressor(nTrialBack+1:end,:);
missFlag = (targVar ==0); targVar(missFlag) = []; regressor(missFlag,:) = [];
regressor = zscore(regressor,0,1);
targVar(targVar==1) = 0; targVar(targVar==-1) = 1;
options = statset('MaxIter',500);
mdl = fitglm(regressor,targVar,'Distribution','binomial','Link','logit','Options',options);
%targVar(targVar==-1) = 2;
%[B, dev, stats] = mnrfit(regressor,targVar);
end

