%% LOAD ANIMALS
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);
reinf = mouseMega.loadReinf;

%% PLOT 1 -- Naive vs. Expert acc and RT
% select 100 trials for acc and RT comparison

tempBin = 200; start = 200;
preLearningTrial = start:start+tempBin-1;
[~,accPre,dpPre] = cellfun(@(x)(fn_getAccBias(x.behav.stimulus(preLearningTrial),x.behav.responseType(preLearningTrial)==1)),...
    mouseMega.mouseCell);

[~,endIdx] = max(reinf.acc,[],1);
postLearningTrial = {};
for i = 1:mouseMega.nMouse
    if endIdx(i) + tempBin/2 > mouseMega.nTrials(i)
        postLearningTrial{i} = (mouseMega.nTrials(i)-tempBin+1): mouseMega.nTrials(i);
    else
        postLearningTrial{i} = endIdx(i) - tempBin/2 + 1 : endIdx(i) + tempBin/2;
    end  
end

fn_plotWheelByStim(mouseMega,10,preLearningTrial,'noColor');
tempIdx = 2200;
fn_plotWheelByStim(mouseMega,10,tempIdx:tempIdx+200,'noColor');
