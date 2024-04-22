%% LOAD DATA
clear; 
mice ={'zz107','zz109','zz111','zz112','zz113','zz115'};
opsParam.biasBlockType = 'threshold';
allMouse = fn_getObjPT_bin30(mice,opsParam); 
mouseMega = wheel2AFCmega(allMouse);
probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
attributesCell= {};animalIdx= {};
for i = 1:mouseMega.nMouse
    attributesCell{i} = allMouse{i}.behav(:,28:33);
    animalIdx{i} = ones(size(attributesCell{i},1),1)*i;
end
tempAttributes = table2array(fn_cell2mat(attributesCell,1));
animalIdx = fn_cell2mat(animalIdx,1);
tempNan = find(isnan(tempAttributes(:,1)));
tempAttributes(tempNan,:) = nanmean(tempAttributes,1);

pythonPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\wheelData\';

save([pythonPath 'wheelAttributes_allMouse_032923.mat'],'tempAttributes','animalIdx','tempNan');