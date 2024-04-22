clear;
mouse = 'zz115';
dataPath = 'C:\Users\zzhu34\Documents\tempData\octoData\trialData\';
targetPath = 'C:\Users\zzhu34\Documents\tempData\octoData\psyTrackData\';
mkdir(targetPath);
load([dataPath filesep mouse '.mat']);
trialData = fn_readStructByFieldKey(trialData,'trainingType',{'puretone'}); 
trialData = fn_selectAnimalDay(trialData,mouse);
nPrev = 5;
nSplit = 5;
removeMissTrial = true;
%% GET TRIAL VARIABLES
stimulus = fn_cell2mat(trialData.stimulus,1); % stimulus, 1 or 2, need to change to -1 and 1
action = fn_cell2mat(trialData.action,1); % action, 1 or 2, need to change to -1 and 1 for action history, 
correct = fn_cell2mat(trialData.responseType,1);% correct or not, 1 or 2, need to chagne to 1 or 0

%% take out miss trials
if removeMissTrial
    missFlag = (action == 0); missIdx = find(missFlag); 
    if mod(sum(~missFlag),nSplit) ~= 0 
        reminder =  mod(sum(~missFlag),nSplit); 
        hitIdx = find(missFlag==0); missFlag(hitIdx(end-(reminder-1):end)) = 1;
    end

    stimulus(missFlag) = []; action(missFlag) = []; correct(missFlag) = [];
else
    missFlag = (action == 0); action(missFlag) = nan; correct(missFlag) = 2; 
end
%% GET DAYLENGTH
allDay_keepMiss = [];
dayLength = [];
tempCountDay = 0;
for i = 1:length(trialData.responseType)
    temp = trialData.responseType{i};  
    if ~isempty(temp)>0; tempCountDay = tempCountDay + 1; end
    allDay_keepMiss = cat(2,allDay_keepMiss,ones(1,length(temp))*tempCountDay);
    temp(temp==0) = []; 
    dayLength(i) = length(temp); 
   
end
dayLength(dayLength == 0 ) = [];


%%
y = action; % Choice, 1 or 2, target variable
answer = stimulus; % Correct Answer (same with stimulus), 1 or 2
correct(correct==2) = 0; % Reward, correct or not. 1 or 0 (changed from 1 and 2)
stimulus = stimulus*2 - 3; % Stimulus, Change from 1 or 2 to -1 or 1;
action = action*2-3; % Choice, Change from 1 or 2 to -1 or 1;
actionXreward = action .* correct; % ChoiceXreward, only positive reward is predictive

stimH = getTrialHistory(stimulus, nPrev); % Stimulus History, i.e. choice X reward
actionH = getTrialHistory(action, nPrev);
actionXrewardH = getTrialHistory(actionXreward, nPrev); % Choice History

sameStimFlag = double(stimH == repmat(stimulus,[1 nPrev])); 
actionXrewardXcurrstimH = actionXrewardH.* sameStimFlag; % PrevStim X Curr Stim (only same stimulus) X PrevChoice X reward

if ~removeMissTrial
    y(missFlag) = []; action(missFlag) = []; stimulus(missFlag) = []; answer(missFlag) = []; correct(missFlag) = [];
    stimH(missFlag,:) = []; actionH(missFlag,:) = []; actionXrewardH(missFlag,:) = []; actionXrewardXcurrstimH(missFlag,:) = [];

    consecutiveMissFlag = sum(isnan(actionH),2)==nPrev;
    
    allDay_keepMiss(missFlag) = [];
    
    y(consecutiveMissFlag) = []; action(consecutiveMissFlag) = []; stimulus(consecutiveMissFlag) = []; answer(consecutiveMissFlag) = []; correct(consecutiveMissFlag) = [];
    stimH(consecutiveMissFlag,:) = []; actionH(consecutiveMissFlag,:) = []; actionXrewardH(consecutiveMissFlag,:) = []; actionXrewardXcurrstimH(consecutiveMissFlag,:) = [];
    
    allDay_keepMiss(consecutiveMissFlag) = [];
end




%% SAVE DATA
if removeMissTrial
save([targetPath filesep mouse '.mat'],'y','answer','correct','stimulus','stimH',...
    'actionH','actionXrewardH','actionXrewardXcurrstimH','dayLength');
else
    dayLength = [];
    for i = 1:max(allDay_keepMiss)
        dayLength(i) = sum(allDay_keepMiss==i);
    end
save([targetPath filesep mouse '_missNan_nPrev' int2str(nPrev) '.mat'],'y','answer','correct','stimulus','stimH',...
    'actionH','actionXrewardH','actionXrewardXcurrstimH','dayLength');
end
%% FUNCTIONS
function trialHistory = getTrialHistory(currTrial, nPrev)

trialHistory = zeros(length(currTrial),nPrev);
for i = 1:nPrev
    trialHistory(:,i) = circshift(currTrial,i);
    trialHistory(1:i,i) = 0; % was nan before for keep miss purpose
end

end