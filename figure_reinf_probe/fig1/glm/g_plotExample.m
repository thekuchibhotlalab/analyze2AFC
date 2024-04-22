clear; 
mouse  = 'zz069';task = 'puretone'; 
allMouse = fn_getObjPT_bin40({mouse}); allMouse = allMouse{1};

%%
tempBin = 100; 

[tempBias, ~] = fn_getAccBiasSmooth(allMouse.behav.stimulus,allMouse.behav.responseType,tempBin);
[~, tempAcc] = fn_getAccBiasSmooth(allMouse.behav.stimulus,allMouse.behav.responseType,tempBin);

figure; subplot(2,1,1);hold on; plot(tempAcc,'Color',[0.71 0.71 0.47],'LineWidth',1.5); plot(allMouse.behav.modelPredAcc,'-','Color',[0.2 0.2 0.2],'LineWidth',1.5); 
ylim([0.25 1.0]); xlim([0 3000]); yticks([0.5 1])
subplot(2,1,2);hold on; plot(tempBias,'Color',[0.71 0.71 0.47],'LineWidth',1.5);plot(allMouse.behav.modelPredBias,'-','Color',[0.2 0.2 0.2],'LineWidth',1.5); 
ylim([-0.7 0.7]); xlim([0 3000]);

tempCorrB = corr(tempBias, allMouse.behav.modelPredBias,'rows','complete');
tempCorrC = corr(tempAcc, allMouse.behav.modelPredAcc,'rows','complete');
[~,tempIdxB] = max(tempCorrB);[~,tempIdxC] = max(tempCorrC);
disp([tempCorrC.^2 tempCorrB.^2])

%% evaluate learning curve correlation of all animals
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
allMouse = fn_getObjPT_bin40(mice);  mouseMega = wheel2AFCmega(allMouse);
tempRsq = zeros(length(allMouse),2);
for i = 1:length(allMouse)

    tempMouse = allMouse{i};
    tempBin = 100; 
    
    [tempBias, tempAcc] = fn_getAccBiasSmooth(tempMouse.behav.stimulus,tempMouse.behav.responseType,tempBin);
    tempCorrB = corr(tempBias, tempMouse.behav.modelPredBias,'rows','complete');
    tempCorrC = corr(tempAcc, tempMouse.behav.modelPredAcc,'rows','complete');
    [~,tempIdxB] = max(tempCorrB);[~,tempIdxC] = max(tempCorrC);
    tempRsq(i,:) = [tempCorrC.^2 tempCorrB.^2];
end 