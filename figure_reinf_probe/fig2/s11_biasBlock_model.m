%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT_bin30(mice);
learningCurveBin = 30;

%%
for i = 1:length(allMouse)
    allMouse{i} = getModelBias (allMouse{i});
end
mouseMega = wheel2AFCmega(allMouse);

plotBlockTransition(mouseMega,'transPlotXlim',50);
%% FUNCTIONS 
function inObj = getModelBias (inObj)

[onL, offL, idxL,~,maxL] = fn_getBlockOnOff(inObj.behav.modelBias>0,inObj.behav.modelBias);
[onR, offR, idxR,~,maxR] = fn_getBlockOnOff(inObj.behav.modelBias<0,inObj.behav.modelBias);

selFlagL = maxL>0.2; selFlagR = maxR>0.2; 
biasBlock.blockL.start = onL(selFlagL); biasBlock.blockL.end = offL(selFlagL);
idxL = idxL(selFlagL); biasBlock.blockL.len = cellfun(@length,idxL,'UniformOutput',true);
biasBlock.blockR.start = onR(selFlagR); biasBlock.blockR.end = offR(selFlagR);
idxR = idxR(selFlagR); biasBlock.blockR.len = cellfun(@length,idxR,'UniformOutput',true);

% calculate the stateFlag for the 1st time
biasBlock.stateFlag = zeros(1,size(inObj.behav,1)); 
for i = 1:length(idxL); biasBlock.stateFlag(idxL{i}) = 1; end
for i = 1:length(idxR); biasBlock.stateFlag(idxR{i}) = -1; end

[biasBlock.blockU.start, biasBlock.blockU.end, idxU,~] = fn_getBlockOnOff(biasBlock.stateFlag==0,biasBlock.stateFlag);
biasBlock.blockU.len =  cellfun(@length,idxU,'UniformOutput',true);
% GET BIAS TRANSITION IDENTITY
[biasBlock.trans, biasBlock.transID] = fn_getTransition(biasBlock);
inObj.biasBlock = biasBlock;

end 