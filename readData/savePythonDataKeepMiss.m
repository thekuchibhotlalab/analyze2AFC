clear;
mouse = 'zz054';
dataPath = 'C:\Users\zzhu34\Documents\tempData\octoData\trialData\';
targetPath = 'C:\Users\zzhu34\Documents\tempData\octoData\psyTrackData\trialData';
mkdir(targetPath);
load([dataPath filesep mouse '.mat']);
trialData = fn_readStructByFieldKey(trialData,'trainingType',{'puretone'});

switch mouse
    case 'zz054'
        maxDay = 20210610;
    case 'zz062'
        maxDay = 20210621;
    case 'zz063'
        maxDay = 20210611;
    
    case 'zz066'
        maxDay = 20210619;
    case 'zz067'
        maxDay = 20210706;
    case 'zz068'
        maxDay = 20210622;
    case 'zz069'
        maxDay = 20210629;
        
    case 'zz111'
        maxDay = 20220312;
    case 'zz112'
        maxDay = 20220312;
    case 'zz113'
        maxDay = 20220312;
    case 'zz115'
        maxDay = 20220312;
            
end

tempDate = cellfun(@str2double,trialData.date);
trialData = fn_selectAnimalDay(trialData,mouse);
%trialData = fn_readStructByFlag(trialData,tempDate<=maxDay);

nPrev = 5;
nSplit = 5;
removeMissTrial = false;
%% GET TRIAL VARIABLES
stimulus = fn_cell2mat(trialData.stimulus,1); % stimulus, 1 or 2, need to change to -1 and 1
action = fn_cell2mat(trialData.action,1); % action, 1 or 2, need to change to -1 and 1 for action history, 
correct = fn_cell2mat(trialData.responseType,1);% correct or not, 1 or 2, need to chagne to 1 or 0

%% GET DAYLENGTH
dayIdx = [];
tempDayLen = [];
tempCountDay = 0;
for i = 1:length(trialData.responseType)
    temp = trialData.responseType{i};  
    if ~isempty(temp)>0; tempCountDay = tempCountDay + 1; tempDayLen(tempCountDay) = length(temp); end
    dayIdx = cat(2,dayIdx,ones(1,length(temp))*tempCountDay); 
end


%%
y = action; % Choice, 1 or 2, target variable
answer = stimulus; % Correct Answer (same with stimulus), 1 or 2
correct(correct==2) = 0; % Reward, correct or not. 1 or 0 (changed from 1 and 2)
stimulus = stimulus*2 - 3; % Stimulus, Change from 1 or 2 to -1 or 1;
action = action*2-3; action(action==-3) = 0;% Choice, Change from 1 or 2 to -1 or 1;
actionXposReward = action .* correct; % ChoiceXreward, only positive reward is predictive
actionXnegReward = action .* (1-correct); % ChoiceXreward, only positive reward is predictive

stimH = getTrialHistory(stimulus, nPrev,tempDayLen); % Stimulus History, i.e. choice X reward
actionH = getTrialHistory(action, nPrev,tempDayLen);
actionXposRewardH = getTrialHistory(actionXposReward, nPrev,tempDayLen); % Choice History
actionXnegRewardH = getTrialHistory(actionXnegReward, nPrev,tempDayLen); % Choice History

sameStimFlag = double(stimH == repmat(stimulus,[1 nPrev])); 
%actionXrewardXcurrstimH = actionXrewardH.* sameStimFlag; % PrevStim X Curr Stim (only same stimulus) X PrevChoice X reward

if ~removeMissTrial
    missFlag = (y==0); 
    %consecutiveMissFlag = sum(isnan(actionH),2)==nPrev;
    
    removeFlag = missFlag;
    keepIdx = find(~removeFlag); 
    removeFlag(keepIdx(end-mod(length(keepIdx),nSplit)+1:end)) = true; 
    
    y(removeFlag) = []; action(removeFlag) = []; stimulus(removeFlag) = []; answer(removeFlag) = []; correct(removeFlag) = [];
    stimH(removeFlag,:) = []; actionH(removeFlag,:) = []; actionXposRewardH(removeFlag,:) = []; actionXnegRewardH(removeFlag,:) = []; %actionXrewardXcurrstimH(missFlag,:) = [];
    
    dayIdx(removeFlag) = [];
end

%% SAVE DATA
for i = 1:max(dayIdx)
    dayLength(i) = sum(dayIdx==i);
end

save([targetPath filesep mouse '_nPrev' int2str(nPrev) '.mat'],'y','answer','correct','stimulus','stimH',...
    'actionH','actionXposRewardH','actionXnegRewardH','dayLength','dayIdx','removeFlag');

%% FUNCTIONS
function trialHistory = getTrialHistory(currTrial, nPrev,tempDayLen)
tempDayLen = [0 cumsum(tempDayLen)];
trialHistory = zeros(length(currTrial),nPrev);
for i = 1:nPrev; trialHistory(:,i) = circshift(currTrial,i); end


for j = 1:length(tempDayLen)-1
    for i = 1:nPrev
        trialHistory(tempDayLen(j)+1:tempDayLen(j)+i,i) = 0;
    end
end

end