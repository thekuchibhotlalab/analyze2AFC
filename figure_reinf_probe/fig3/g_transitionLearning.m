%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT_bin30(mice);
learningCurveBin = 30;
mouseMega = wheel2AFCmega(allMouse);

%%
plotBlockTransitionFreq(mouseMega)