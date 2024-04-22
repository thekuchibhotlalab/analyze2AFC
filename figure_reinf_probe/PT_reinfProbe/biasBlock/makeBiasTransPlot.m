%%
clear;
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);
allMouse = fn_getObjPT_bin30(mice);
mouseMega = wheel2AFCmega(allMouse);

%%
plotBlockTransition(mouseMega)
%%
plotBlockTransitionFreq(mouseMega)

%%
mouseRand = rand(1,1000000); mouseRand(mouseRand<0.5) = -1; mouseRand(mouseRand>0.5) = 1;
biasRand = smoothdata(mouseRand,'movmean',30);
[~, ~, ~,bareaL] = fn_getBlockOnOff(biasRand>0,biasRand);
[~, ~, ~,bareaR] = fn_getBlockOnOff(biasRand<0,biasRand);
barea = [bareaL bareaR]; areaThre = prctile(barea,95);

behavArea = [];
for i = 1:mouseMega.nMouse
    tempBias = mouseMega.mouseCell{i}.behav.bias; 
    [~, ~, ~,areaL] = fn_getBlockOnOff(tempBias>0,tempBias);
    [~, ~, ~,areaR] = fn_getBlockOnOff(tempBias<0,tempBias);
    behavArea = [behavArea areaL areaR];
end

startVal = areaThre/256;
figure; hold on;
plot(log([areaThre areaThre]),[0 0.3],'Color',[0.6 0.6 0.6],'LineWidth',1)
h1 = fn_plotHistLine(barea,'histCountArgIn',{startVal*2.^([0 1 2:0.5:13]),'Normalization','probability'},'xaxislog',true);
set(h1,'Color',[0.6 0.6 0.6]); set(h1,'LineWidth',2);
h2 = fn_plotHistLine(behavArea,'histCountArgIn',{startVal*2.^([0 1 2:0.5:13]),'Normalization','probability'},'xaxislog',true);
set(h2,'Color',fn_wheelColorsPT('bias')); set(h2,'LineWidth',2);
ylim([0 0.3])
xticks(log([0.1 1 10 100])); xticklabels([0.1 1 10 100]); xlim([log(startVal)-0.3 log(startVal*2^13)+0.3]);
ylabel('Probability'); xlabel('AUC of each bias block')


startVal = areaThre;
figure; hold on;
plot(log([areaThre areaThre]),[0 0.3],'Color',[0.6 0.6 0.6],'LineWidth',1)
h1 = fn_plotHistLine(barea,'histCountArgIn',{startVal*2.^(0:0.25:6),'Normalization','probability'},'xaxislog',true);
set(h1,'Color',[0.6 0.6 0.6]); set(h1,'LineWidth',2);
h2 = fn_plotHistLine(behavArea,'histCountArgIn',{startVal*2.^(0:0.25:6),'Normalization','probability'},'xaxislog',true);
set(h2,'Color',fn_wheelColorsPT('bias')); set(h2,'LineWidth',2);
ylim([0 0.02])
xticks(log([10 40 160 640])); xticklabels([10 40 160 640]); xlim([log(startVal)-0.3 log(startVal*2^6)+0.3]);
ylabel('Probability'); xlabel('AUC of each bias block')
