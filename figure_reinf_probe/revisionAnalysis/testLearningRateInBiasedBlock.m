clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT_bin40(mice);
learningCurveBin = 40;
%%
accGainL = {}; accGainR = {}; accGainU = {}; accGainUL = {}; accGainUR = {};
aRU = {}; aRL = {}; aRR = {};
trialBefAft = 40;
for i = 1:13
    learningStart = find(allMouse{i}.behav.modelAcc>0.5,1,'first'); 
    learningEnd = find(allMouse{i}.behav.modelAcc>0.8,1,'first'); 
    for j = 1:length(allMouse{i}.biasBlock.blockL.start)
        
        if  allMouse{i}.biasBlock.blockL.start(j) >= learningStart && allMouse{i}.biasBlock.blockL.start(j) <= learningEnd
            tempStart =  allMouse{i}.biasBlock.blockL.start(j);
            tempEnd = allMouse{i}.biasBlock.blockL.end(j);
            accGainL{i}(j) = (allMouse{i}.behav.modelAcc(tempEnd) - allMouse{i}.behav.modelAcc(tempStart))...
                / (tempEnd-tempStart+1);
            try 
                temp1 = (allMouse{i}.behav.modelAcc(tempEnd+trialBefAft) - allMouse{i}.behav.modelAcc(tempEnd))/ (trialBefAft+1);
            catch
                temp1 = nan;
            end 
                try
                temp2 = (allMouse{i}.behav.modelAcc(tempStart) - allMouse{i}.behav.modelAcc(tempStart-trialBefAft))/ (trialBefAft+1);
            catch 
                temp2 = nan; 
            end 
            accGainUL{i}(j) = nanmean([temp1 temp2]);

        else 
            accGainL{i}(j) = nan; 
        end
    end

    for j = 1:length(allMouse{i}.biasBlock.blockR.start)
        if allMouse{i}.biasBlock.blockR.start(j) >= learningStart && allMouse{i}.biasBlock.blockR.start(j) <= learningEnd
            tempStart =  allMouse{i}.biasBlock.blockR.start(j);
            tempEnd = allMouse{i}.biasBlock.blockR.end(j);
            accGainR{i}(j) = (allMouse{i}.behav.modelAcc(tempEnd) - allMouse{i}.behav.modelAcc(tempStart))...
                / (tempEnd-tempStart+1);
            try
                temp1 = (allMouse{i}.behav.modelAcc(tempEnd+trialBefAft) - allMouse{i}.behav.modelAcc(tempEnd))/ (trialBefAft+1);
            catch
                temp1 = nan; 
            end
            try
                temp2 = (allMouse{i}.behav.modelAcc(tempStart) - allMouse{i}.behav.modelAcc(tempStart-trialBefAft))/ (trialBefAft+1);
            catch 
                temp2 = nan; 
            end
            accGainUR{i}(j) = nanmean([temp1 temp2]);
        else
            accGainR{i}(j) = nan; 
        end 
    end

    for j = 1:length(allMouse{i}.biasBlock.blockU.start)
        if allMouse{i}.biasBlock.blockU.start(j) >= learningStart && allMouse{i}.biasBlock.blockU.start(j) <= learningEnd
            accGainU{i}(j) = allMouse{i}.behav.modelAcc(allMouse{i}.biasBlock.blockU.end(j)) - ...
                allMouse{i}.behav.modelAcc(allMouse{i}.biasBlock.blockU.start(j));
            accGainU{i}(j) = accGainU{i}(j) / allMouse{i}.biasBlock.blockU.len(j);
        else
            accGainU{i}(j) = nan; 
        end
    end
end

accGainU = cellfun(@nanmean,accGainU);
accGainBias = cellfun(@(x,y)(nanmean([x y])),accGainL,accGainR);
accGainBiasU = cellfun(@(x,y)(nanmean([x y])),accGainUL,accGainUR);

scaleFactor = 10000;
figure; 
[p] = fn_plotComparison({accGainBias*scaleFactor, accGainU*scaleFactor},'paired',true);
xticks([1 2]); xticklabels({'Biased epoch','Unbiased epoch'})
ylabel('Accuracy gain (%) per 100 trials')
figure; 
[p] = fn_plotComparison({accGainBias*scaleFactor, accGainBiasU*scaleFactor},'paired',true);
xticks([1 2]); xticklabels({'Biased epoch','Unbiased epoch'});
ylabel('Accuracy gain (%) per 100 trials')