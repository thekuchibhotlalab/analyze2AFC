clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT_bin40(mice);
learningCurveBin = 40;
%% 
mouseRand = rand(1,10000000); mouseRand(mouseRand<0.5) = -1; mouseRand(mouseRand>0.5) = 1;
biasRand = smoothdata(mouseRand,'movmean',learningCurveBin);
[~, ~, ~,bareaL] = fn_getBlockOnOff(biasRand>0,biasRand);
[~, ~, ~,bareaR] = fn_getBlockOnOff(biasRand<0,biasRand);
barea = [bareaL bareaR]; areaThre = prctile(barea,95);

biasAUC = {};
for i = 1:length(allMouse)
    biasAUC{i} = [];
    for j = 1:length(allMouse{i}.biasBlock.blockL.start)
        tempIdx = allMouse{i}.biasBlock.blockL.start(j):allMouse{i}.biasBlock.blockL.end(j);
        biasAUC{i} = [biasAUC{i}, nansum(abs(allMouse{i}.behav.bias(tempIdx)))];
    end 

    for j = 1:length(allMouse{i}.biasBlock.blockR.start)
        tempIdx = allMouse{i}.biasBlock.blockR.start(j):allMouse{i}.biasBlock.blockR.end(j);
        biasAUC{i} = [biasAUC{i}, nansum(abs(allMouse{i}.behav.bias(tempIdx)))];
    end 

end
biasAUC = fn_cell2mat(biasAUC,2);
binExp = areaThre * 2.^(-9:1:8);
figure; hold on; 

fn_plotHistLine(barea,'histCountArgIn',{binExp,'Normalization','probability'},...
    'plotArgIn',{'LineWidth',2,'Color',[0.8 0.8 0.8]});

fn_plotHistLine(biasAUC,'histCountArgIn',{binExp,'Normalization','probability'},...
    'plotArgIn',{'LineWidth',2,'Color',fn_wheelColorsPT('bias')});
plot([areaThre areaThre],[0 0.5],'Color',[0.2 0.2 0.2, 0.6],'LineWidth',2)

xlim([binExp(1) binExp(end)]); ylim([0 0.5]); xticks([0.1 1 10 100 1000])
set(gca,'XScale','log')