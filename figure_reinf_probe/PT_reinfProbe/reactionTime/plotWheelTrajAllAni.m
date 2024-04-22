%% PART 1 --CLUSTERING

%% PART 1.1 -- GET PARAMETERS FOR CLUSTERING
clear; 
mice ={'zz107','zz109','zz111','zz112','zz113','zz115'};
allMouse = fn_getObjPT_bin30(mice);mouseMega = wheel2AFCmega(allMouse);
%%
saveWheelTrajAttr(mouseMega, mice);

%save the data of all animals
tempAttributes = []; animalIdx = []; 
for i = 1:6; a = load(['wheelAttributes_' mice{i} '.mat']);
animalIdx( length(animalIdx) + (1:length(a.tempAttributes)) )  = i;
tempAttributes = cat(1,tempAttributes,a.tempAttributes);
end
save('wheelAttributes_allMouse.mat','tempAttributes','animalIdx');
%% PART 2.1 -- LOAD DATA FROM PYTHON
clear; 
pythonPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\wheelData\';

load([pythonPath filesep 'wheelAttributes_allMouse.mat' ]);
load([pythonPath filesep 'wheel_cluster_allMouse.mat' ]); 
load([pythonPath filesep 'wheel_cluster_allMouse_tsne.mat' ]);

mice ={'zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice); opsParam.biasBlockType = 'threshold';
allMouse = fn_getObjPT_bin30(mice,opsParam);
mouseMega = wheel2AFCmega(allMouse);

label4 = double(label4 + 1); label5 = double(label5 + 1); 
attributesOrig = tempAttributes;

%% Plot 
allRT = fn_cell2mat(cellfun(@(x)(x.behav.reactionTime),allMouse,'UniformOutput',false),1);


figure; subplot(2,2,1);hold on;
for i = 1:4
    tempFlag = label4==i;
    scatter3(allRT(tempFlag), attributesOrig(tempFlag,1),attributesOrig(tempFlag,6),10,matlabColors(i),'filled','MarkerFaceAlpha',0.3)
end
xlabel('reaction time'); ylabel('onset time'); zlabel('total distance'); 



figure; subplot(1,3,1);hold on;
for i = 1:4
    tempFlag = label4==i;
    scatter(attributesOrig(tempFlag,1),attributesOrig(tempFlag,6),10,matlabColors(i),'filled','MarkerFaceAlpha',0.3)
end
xlabel('onset time'); ylabel('total distance')

subplot(1,3,2);hold on;
for i = 1:4
    tempFlag = label4==i;
    cdfplot(attributesOrig(tempFlag,1))
end
xlabel('onset time'); ylabel('cumulative distribution')

subplot(1,3,3);hold on;
for i = 1:4
    tempFlag = label4==i;
    cdfplot(attributesOrig(tempFlag,6))
end
xlabel('onset time'); ylabel('cumulative distribution')
%% PART 2.2 -- PLOT CLUSTERING SUMMARY AND CLUSTER SELECTION

% PLOT CLSUTERING RESULTS FOR ALL ANIMALS

[wheelDownSample, RT,selFlag] = loadWheelClusterData(mouseMega, []);

clusterProportion = plotCluster(wheelDownSample,RT,attributesScaled(selFlag,:),tempAttributes(selFlag,:), animalIdx,mice,label4(selFlag),4);
%clusterProportion = plotCluster(wheelDownSample_all,RT_all,attribute_tsne,tempAttributes, animalIdx,mice,label5,5);
figure; hold on; plot(cluster_metric,'Color',[0.2 0.2 0.2],'LineWidth',2); 
scatter(1:length(cluster_metric),cluster_metric,30,[0.2 0.2 0.2],'filled')
ylabel('Sum of squared error (WSS)'); xlabel('number of clusters'); 

figure; subplot(2,1,1); hold on
fastProp = clusterProportion{1} + clusterProportion{2};
f_errorbar = fn_plotFillErrorbar(1:size(fastProp,2),nanmean(fastProp),fn_nansem(fastProp),...
        matlabColors(1),'faceAlpha',0.3,'LineStyle','none');
plot(1:size(fastProp,2),nanmean(fastProp),'Color',matlabColors(1),'LineWidth',2);
xlim([0 2500]); ylim([0.5 1]);yticks(0.5:0.25:1)
subplot(2,1,2); hold on
slowProp = clusterProportion{3} + clusterProportion{4};
f_errorbar = fn_plotFillErrorbar(1:size(slowProp,2),nanmean(slowProp),fn_nansem(slowProp),...
        matlabColors(1),'faceAlpha',0.3,'LineStyle','none');
plot(1:size(slowProp,2),nanmean(slowProp),'Color',matlabColors(1),'LineWidth',2);
xlim([0 2500]); ylim([0 0.5]); yticks(0:0.25:0.5)
%% PART 2.3 -- SEPARATE CLUSTER DATA FOR ALL ANIMALS
[wheelDownSample, RT,selFlag] = loadWheelClusterData(mouseMega, []);
bias = fn_cell2mat(mouseMega.getProp('behav','field','bias','matFlag',false),1);
plotBiasRT(RT,bias,label4);


%% PART 2.4 -- SEPARATE CLUSTER DATA FOR EACH ANIMAL
selMouse = 1; selDay = 2:4;
bias = fn_cell2mat(mouseMega.getProp('behav','field','bias','matFlag',false,'idx',selMouse),1);
day = fn_cell2mat(mouseMega.getProp('behav','field','day','matFlag',false,'idx',selMouse),1);
[wheelDownSample, RT,selFlag] = loadWheelClusterData(mouseMega, selMouse);
dayFlag = sum(day==selDay,2) > 0;
bias = bias(dayFlag);wheelDownSample = wheelDownSample(dayFlag,:); RT = RT(dayFlag);
label = label4(selFlag); label = label(dayFlag);

plotBiasRT(RT,bias,label);

%% PLOT REACTION TIME 

[wheelDownSample, RT,selFlag] = loadWheelClusterData(mouseMega, []);
action = fn_cell2mat(mouseMega.getProp('behav','field','action','matFlag',false),1);
plotActionRT(RT,action,label4);

%%
selMouse = 6; 
action = fn_cell2mat(mouseMega.getProp('behav','field','action','matFlag',false,'idx',selMouse),1);
day = fn_cell2mat(mouseMega.getProp('behav','field','day','matFlag',false,'idx',selMouse),1); disp([max(day)])
selDay = (max(day)-5):max(day);
[wheelDownSample, RT,selFlag] = loadWheelClusterData(mouseMega, selMouse);
dayFlag = sum(day==selDay,2) > 0;
action = action(dayFlag);wheelDownSample = wheelDownSample(dayFlag,:); RT = RT(dayFlag);
label = label4(selFlag); label = label(dayFlag);

plotActionRT(RT,action,label)



%% PART 2.6 -- PLOT LM RESULTS FOR EACH ANIMAL
LM_rsq = table2array(tblLM_rsq); figure; bar(mean(LM_rsq,2));
xticklabels(tblLM_rsq.Properties.RowNames);

biasP = fn_cell2mat(cellfun(@(x)(x.bias<0.05),tblLM_coeffP,'UniformOutput',false),2); 
figure; bar(mean(biasP,2));xticklabels(tblLM_rsq.Properties.RowNames);

cbiasP = fn_cell2mat(cellfun(@(x)(x.choiceXbias<0.05),tblLM_coeffP,'UniformOutput',false),2); 
figure; bar(mean(cbiasP,2));xticklabels(tblLM_rsq.Properties.RowNames);

rewardP = fn_cell2mat(cellfun(@(x)(x.rewardH1<0.05),tblLM_coeffP,'UniformOutput',false),2); 
figure; bar(mean(rewardP,2));xticklabels(tblLM_rsq.Properties.RowNames);

ctxtP = fn_cell2mat(cellfun(@(x)(x.context<0.05),tblLM_coeffP,'UniformOutput',false),2); 
figure; bar(mean(ctxtP,2));xticklabels(tblLM_rsq.Properties.RowNames);

%% PART 3.0 -- COSNTRUCT LINEAR REGRESSION MODEL FOR ALL ANIMALS COMBINED
% THIS DOES NOT WORK SINCE DIFFERENT ANIMALS ARE BIASED TO DIFFERENT SIDES
%tblLM_rsq = table(); tblLM_coeff = table(); tblLM_coeffP = table();
%[tblLM_rsq,tblLM_coeff, tblLM_coeffP] = fitLM1(mouseMega,[],label4, attributesScaled);


%% PART 3.0 -- COSNTRUCT LINEAR REGRESSION MODEL FOR EACH ANIMAL, ALL DAYS
selMouse = num2cell(1:6);
[tblCell,RTCell,labelCell,attributeCell,~] = cellfun(@(x)(makeRegressorTable(mouseMega,x,label4,attributesScaled)),selMouse,'UniformOutput',false);
[tblLM_rsq,tblLM_coeff, tblLM_coeffP, tblLM_uniContr,pred] = cellfun(@(a,b,c,d)(fitLM_allVar(a,b,c,d,true)), tblCell...
    ,RTCell,labelCell,attributeCell, 'UniformOutput', false);
temp = tblLM_rsq;tblLM_rsq = table();
for i = 1:mouseMega.nMouse; tblLM_rsq(:,i) = temp{i}; end 


%% PART 3.1 -- VISUALIZATION OF 3.0
varNames = tblLM_uniContr{1}.Properties.VariableNames;
uniVarExp = table();
for i = 1:length(varNames)
    tempMat = fn_cell2mat(cellfun(@(x)(x.(varNames{i})),tblLM_uniContr,'UniformOutput',false),2);
    uniVarExp.(varNames{i}) = mean(tempMat ./ table2array(tblLM_rsq),2)*100;
end
uniVarExp.Properties.RowNames = tblLM_uniContr{1}.Properties.RowNames;

figure; subplot(1,2,1); imagesc(table2array(tblLM_rsq)); colorbar; caxis([0 0.2]) 
xticks(1:size(uniVarExp,2)); xticklabels(mouseMega.mouseName)
yticks(1:size(uniVarExp,1)); yticklabels(tblLM_rsq.Properties.RowNames)
subplot(1,2,2);imagesc(table2array(uniVarExp)); colorbar; caxis([0 20]); 
xticks(1:size(uniVarExp,2)); xticklabels(uniVarExp.Properties.VariableNames)
yticks(1:size(uniVarExp,1)); yticklabels(uniVarExp.Properties.RowNames)

%%
selMouse = 3;
[~, RT,selFlag] = loadWheelClusterData(mouseMega, selMouse);
templabel = label4(selFlag); RT_fast = RT(templabel<=2,:); RT_slow = RT(templabel>2,:);
figure; subplot(2,2,1);plot([0 2.5],[0 2.5],'Color',[0.8 0.8 0.8]); hold on; scatter(RT,pred{selMouse}{1},10,'filled'); xlim([0 2.5]); ylim([0 2.5]);
title('Regression of RT'); xlabel('RT'); ylabel('predicted RT')
subplot(2,2,2);plot([0 2.5],[0 2.5],'Color',[0.8 0.8 0.8]); hold on; scatter(log(RT),pred{selMouse}{4},10,'filled'); xlim([0 2.5]); ylim([0 2.5]);
title('Regression of RT'); xlabel('RT'); ylabel('predicted RT')
subplot(2,2,3); plot([0 2.5],[0 2.5],'Color',[0.8 0.8 0.8]); hold on; scatter(RT_fast,pred{selMouse}{2},10,'filled'); xlim([0 1.5]); ylim([0 1.5])
title('Regression of RT fast'); xlabel('RT'); ylabel('predicted RT')
subplot(2,2,4); plot([0 2.5],[0 2.5],'Color',[0.8 0.8 0.8]); hold on; scatter(RT_slow,pred{selMouse}{3},10,'filled'); xlim([0.5 2.5]); ylim([0.5 2.5])
title('Regression of RT slow'); xlabel('RT'); ylabel('predicted RT')
%% PART 3.2 -- REGRESSION OF REINF VS. PROBE
% IDENTIFY THE LEARNING ONSET OF EACH ANIMAL
probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
outCell = mouseMega.objFun('binProbeByTrialFromLearningOnset',{[reinfThre probeThre],probeTrialBin});
learningOnset = cell2mat(outCell{1,3}); learningStart = learningOnset - 800; learningEnd = learningOnset + 800;
learningStart(learningStart<=0) = 1; 

tblLM_rsq = {}; tblLM_coeff = {}; tblLM_coeffP = {}; tblLM_uniContr = {}; pred = {}; selRT = {};
probeWheel = {}; reinfWheel = {};

for i = 1:mouseMega.nMouse
    reinfIdx = find(logical(mouseMega.mouseCell{i}.behav.reinfBef) | logical(mouseMega.mouseCell{i}.behav.reinfAft));
    reinfIdx(reinfIdx<learningStart(i)) = []; reinfIdx(reinfIdx>learningEnd(i)) = [];
    probeIdx = find(logical(mouseMega.mouseCell{i}.behav.probe));
    probeIdx(probeIdx<learningStart(i)) = []; probeIdx(probeIdx>learningEnd(i)) = [];

    [tbl_ctxt,tempRT,tempLabel,tempAttributes,~] = makeRegressorTable(mouseMega,i,label4,attributesScaled); 

    tbl_ctxt = tbl_ctxt([probeIdx; reinfIdx],:); tempRT = tempRT([probeIdx; reinfIdx]);
    tempLabel = tempLabel([probeIdx; reinfIdx]); tempAttributes = tempAttributes([probeIdx; reinfIdx],:);

    [tblLM_rsq{i},tblLM_coeff{i}, tblLM_coeffP{i}, tblLM_uniContr{i},pred{i}] = fitLM_allVar(tbl_ctxt,tempRT, tempLabel, tempAttributes);
    selRT{i} = tempRT;
end
temp = tblLM_rsq;tblLM_rsq = table();
for i = 1:mouseMega.nMouse; tblLM_rsq(:,i) = temp{i}; end 

%% PART 3.3 -- VISUALIZATION OF 3.2
varNames = tblLM_uniContr{1}.Properties.VariableNames;
uniVarExp = table();
for i = 1:length(varNames)
    tempMat = fn_cell2mat(cellfun(@(x)(x.(varNames{i})),tblLM_uniContr,'UniformOutput',false),2);
    uniVarExp.(varNames{i}) = mean(tempMat ./ table2array(tblLM_rsq),2)*100;
end
uniVarExp.Properties.RowNames = tblLM_uniContr{1}.Properties.RowNames;

figure; subplot(1,2,1); imagesc(table2array(tblLM_rsq)); colorbar; caxis([0 0.2]) 
xticks(1:size(uniVarExp,2)); xticklabels(mouseMega.mouseName)
yticks(1:size(uniVarExp,1)); yticklabels(tblLM_rsq.Properties.RowNames)
subplot(1,2,2);imagesc(table2array(uniVarExp)); colorbar; caxis([0 20]); 
xticks(1:size(uniVarExp,2)); xticklabels(uniVarExp.Properties.VariableNames)
yticks(1:size(uniVarExp,1)); yticklabels(uniVarExp.Properties.RowNames)

selMouse = 3;
figure; plot([0 2.5],[0 2.5],'Color',[0.8 0.8 0.8]); hold on; scatter(selRT{selMouse},pred{selMouse}{1},10,'filled'); xlim([0 2.5]); ylim([0 2.5]);
title('Regression of RT'); xlabel('RT'); ylabel('predicted RT')

%% PART 3.4 -- COMPARE LINEAR REGRESSION MODELS ACROSS ANIMALS 
selMouse = num2cell(1:6);

[tblCell,RTCell,labelCell,attributeCell,~] = cellfun(@(x)(makeRegressorTable(mouseMega,x,label4,attributesScaled)),selMouse,'UniformOutput',false);
selDay = 4:7; 

RTCell = cellfun(@(x,y)( y((sum(x.day==selDay,2))>0,:)),tblCell,RTCell,'UniformOutput',false);
labelCell = cellfun(@(x,y)( y((sum(x.day==selDay,2))>0,:)),tblCell,labelCell,'UniformOutput',false);
attributeCell = cellfun(@(x,y)( y((sum(x.day==selDay,2))>0,:)),tblCell,attributeCell,'UniformOutput',false);
tblCell = cellfun(@(x,y)( y((sum(x.day==selDay,2))>0,:)),tblCell,tblCell,'UniformOutput',false);

removeVarNames = {'context','rewardH3','rewardH2','day','trials','bias','choiceXbias','actionRate'};
tblLM_rsq = []; AIC = [];
for i = 1:(length(removeVarNames)+1)
    if i~=1
        for j = 1:length(tblCell); tblCell{j}.(removeVarNames{i-1}) = []; end
    end
    [tblLM_rsq_temp,~, ~, ~,AIC_temp, ~] = cellfun(@(a,b,c,d)(fitLM_allVar(a,b,c,d,false)), tblCell...
        ,RTCell,labelCell,attributeCell, 'UniformOutput', false);
    temp = tblLM_rsq_temp;tblLM_rsq_temp = table();
    for k = 1:mouseMega.nMouse; tblLM_rsq_temp(:,k) = temp{k}; end; tblLM_rsq(:,:,i) = table2array(tblLM_rsq_temp);
    temp = AIC_temp;AIC_temp = table();
    for k = 1:mouseMega.nMouse; AIC_temp(:,k) = temp{k}; end; AIC(:,:,i) = table2array(AIC_temp);
end
%% PART 3.5 -- COMPARE LINEAR REGRESSION MODELS ACROSS ANIMALS 
modelNames = {'+context','+rewardH3','+rewardH2','+day','+trials','+bias','+choiceXbias','+actionRate','action-only'};
plotModelComparison(AIC,tblLM_rsq_temp.Properties.RowNames,modelNames,'AIC');
plotModelComparison(tblLM_rsq,tblLM_rsq_temp.Properties.RowNames,modelNames,'Var Exp');

selModel = 6; selModelFit = tblLM_rsq(:,:,selModel);
figure;  fn_plotBar(selModelFit); xticklabels(tblLM_rsq_temp.Properties.RowNames); ylabel('R-squared')

selModel = 4; selModelFit = tblLM_rsq(:,:,selModel);
figure;  fn_plotBar(selModelFit); xticklabels(tblLM_rsq_temp.Properties.RowNames); ylabel('R-squared')
%% PART 5.0 - WHEEL TRAJECTORY REINF VS. PROBE VISUALIZATION

probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
outCell = mouseMega.objFun('binProbeByTrialFromLearningOnset',{[reinfThre probeThre],probeTrialBin});
learningOnset = cell2mat(outCell{1,3}); learningStart = learningOnset - 800; learningEnd = learningOnset + 800;
learningStart(learningStart<=0) = 1; 
probeWheel = {}; reinfWheel = {}; allTbl_ctxt = {};
for i = 1:mouseMega.nMouse
    reinfIdx = find(logical(mouseMega.mouseCell{i}.behav.reinfBef));% | logical(mouseMega.mouseCell{i}.behav.reinfAft));
    reinfIdx(reinfIdx<learningStart(i)) = []; reinfIdx(reinfIdx>learningEnd(i)) = [];
    probeIdx = find(logical(mouseMega.mouseCell{i}.behav.probe));
    probeIdx(probeIdx<learningStart(i)) = []; probeIdx(probeIdx>learningEnd(i)) = [];
    [tbl_ctxt,tempRT,tempLabel,tempAttributes,~] = makeRegressorTable(mouseMega,i,label4,attributesOrig); 
    tbl_ctxt.RT = tempRT; tbl_ctxt.context(tbl_ctxt.context==3) = 2;
    tbl_ctxt.onTime = tempAttributes(:,1); tbl_ctxt.onSpeed = tempAttributes(:,2); tbl_ctxt.onDist = tempAttributes(:,3);
    tbl_ctxt.totalTime = tempAttributes(:,4); tbl_ctxt.totalSpeed = tempAttributes(:,5); tbl_ctxt.totalDist = tempAttributes(:,6);

    tbl_ctxt = tbl_ctxt([probeIdx; reinfIdx],:); %tempRT = tempRT([probeIdx; reinfIdx]); tempAttributes = tempAttributes([probeIdx; reinfIdx],:);
    

    [wheelDownSample, ~,selFlag] = loadWheelClusterData(mouseMega, i);
    templabel = label4(selFlag);  tbl_ctxt.label = templabel([probeIdx; reinfIdx])';
    allTbl_ctxt{i} = tbl_ctxt; 

    temp = wheelDownSample([probeIdx; reinfIdx],:);
    probeWheel{i,1} = wheelDownSample(probeIdx,:); reinfWheel{i,1} = wheelDownSample(reinfIdx,:);
    probeWheel{i,2} = temp(tbl_ctxt.action ==1 & tbl_ctxt.context ==2 & tbl_ctxt.label <=2,:); 
    probeWheel{i,3} = temp(tbl_ctxt.action ==2 & tbl_ctxt.context ==2 & tbl_ctxt.label <=2, :);
    reinfWheel{i,2} = temp(tbl_ctxt.action ==1 & tbl_ctxt.context ==1 & tbl_ctxt.label <=2,:); 
    reinfWheel{i,3} = temp(tbl_ctxt.action ==2 & tbl_ctxt.context ==1 & tbl_ctxt.label <=2,:);
end
%%
selMouse = 1;
figure; subplot(1,2,1); plot(probeWheel{selMouse,1}'); title('probe');hold on; subplot(1,2,2); plot(reinfWheel{selMouse,1}');title('reinf');
figure; subplot(1,2,1); hold on; 
fn_plotFillErrorbar(1:size(probeWheel{selMouse,2},2),nanmean(probeWheel{selMouse,2},1),...
    nanstd(probeWheel{selMouse,2})./sqrt(size(probeWheel{selMouse,2},1)),matlabColors(1,0.2),'LineStyle','none');

fn_plotFillErrorbar(1:size(reinfWheel{selMouse,2},2),nanmean(reinfWheel{selMouse,2},1),...
    nanstd(reinfWheel{selMouse,2})./sqrt(size(reinfWheel{selMouse,2},1)),matlabColors(2,0.2),'LineStyle','none');
plot(nanmean(probeWheel{selMouse,2},1),'Color',matlabColors(1));
plot(nanmean(reinfWheel{selMouse,2},1),'Color',matlabColors(2)); title('Left');

subplot(1,2,2); hold on; 
fn_plotFillErrorbar(1:size(probeWheel{selMouse,3},2),nanmean(probeWheel{selMouse,3},1),...
    nanstd(probeWheel{selMouse,3})./sqrt(size(probeWheel{selMouse,3},1)),matlabColors(1,0.2),'LineStyle','none');
fn_plotFillErrorbar(1:size(reinfWheel{selMouse,3},2),nanmean(reinfWheel{selMouse,3},1),...
    nanstd(reinfWheel{selMouse,3})./sqrt(size(reinfWheel{selMouse,3},1)),matlabColors(2,0.2),'LineStyle','none');
plot(nanmean(probeWheel{selMouse,3},1),'Color',matlabColors(1));
plot(nanmean(reinfWheel{selMouse,3},1),'Color',matlabColors(2)); title('Right');

tempFlag = allTbl_ctxt{selMouse}.context ==1; 
plotActionRT(allTbl_ctxt{selMouse}.RT(tempFlag),allTbl_ctxt{selMouse}.action(tempFlag),allTbl_ctxt{selMouse}.label(tempFlag)',{'L','R'});

tempFlag = allTbl_ctxt{selMouse}.context ==2;
plotActionRT(allTbl_ctxt{selMouse}.RT(tempFlag),allTbl_ctxt{selMouse}.action(tempFlag),allTbl_ctxt{selMouse}.label(tempFlag)',{'L','R'});

tempFlag = allTbl_ctxt{selMouse}.action ==1; 
plotActionRT(allTbl_ctxt{selMouse}.RT(tempFlag),allTbl_ctxt{selMouse}.context(tempFlag),allTbl_ctxt{selMouse}.label(tempFlag)',{'reinf','probe'});

tempFlag = allTbl_ctxt{selMouse}.action ==2;
plotActionRT(allTbl_ctxt{selMouse}.RT(tempFlag),allTbl_ctxt{selMouse}.context(tempFlag),allTbl_ctxt{selMouse}.label(tempFlag)',{'reinf','probe'});

%% PART 5.1 - RT REINF VS. PROBE FAST VS. SLOW CLUSTERS
 [reinfRT,probeRT] = plotReinfProbeDiff(allTbl_ctxt,'RT');
%% PART 5.2 - WHEEL ATTRIBUTES REINF VS. PROBE
plotReinfProbeDiff(allTbl_ctxt,'onTime');
plotReinfProbeDiff(allTbl_ctxt,'onSpeed');
plotReinfProbeDiff(allTbl_ctxt,'onDist');
%%
plotReinfProbeDiff(allTbl_ctxt,'totalTime');
plotReinfProbeDiff(allTbl_ctxt,'totalSpeed');
plotReinfProbeDiff(allTbl_ctxt,'totalDist');

%% PLOT 4.0 -- GET PERIOD OF BIAS VS UNBIASED STATE, for old method of bias block detection
probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
outCell = mouseMega.objFun('binProbeByTrialFromLearningOnset',{[reinfThre probeThre],probeTrialBin});
learningOnset = cell2mat(outCell{1,3}); learningStart = learningOnset - 800; learningEnd = learningOnset + 800;
learningStart(learningStart<=100) = 100; 

nTrialBefAftBlock = 20;
trialNum_bias = cell(mouseMega.nMouse,3); allTbl_bias = {}; allTbl_regression = {}; allWheel = {};
RTcell = {};labelCell = {};attributeCell = {};
for i = 1:mouseMega.nMouse
    tempL = allMouse{i}.biasBlock.blockL.start; tempL_len = allMouse{i}.biasBlock.blockL.len; 
    tempR = allMouse{i}.biasBlock.blockR.start; tempR_len = allMouse{i}.biasBlock.blockR.len; 
    for j = 1:length(tempL)
        temp = min([tempL_len(j) nTrialBefAftBlock]);
        trialNum_bias{i,2} = cat(2,trialNum_bias{i,2},tempL(j)-temp:tempL(j)+temp);  
    end
    for j = 1:length(tempR)
        temp = min([tempR_len(j) nTrialBefAftBlock]);
        trialNum_bias{i,3} = cat(2,trialNum_bias{i,3},tempR(j)-temp:tempR(j)+temp);
    end
    trialNum_bias{i,2} = unique(trialNum_bias{i,2}); trialNum_bias{i,3} = unique(trialNum_bias{i,3});
    trialNum_bias{i,2}( trialNum_bias{i,2}<1) = []; trialNum_bias{i,3}( trialNum_bias{i,3}<1) = [];
    %trialNum_bias{i,2}(trialNum_bias{i,2}<learningStart(i) | trialNum_bias{i,2}>learningEnd(i)) = [];
    %trialNum_bias{i,3}(trialNum_bias{i,3}<learningStart(i) | trialNum_bias{i,3}>learningEnd(i)) = [];
    trialNum_bias{i,1} = cat(2, trialNum_bias{i,2}, trialNum_bias{i,3});


    [tbl_ctxt,tempRT,tempLabel,tempAttributes,~] = makeRegressorTable(mouseMega,i,label4,attributesOrig); 
    tbl_regression = tbl_ctxt;
    tbl_regression.bias = abs(allMouse{i}.biasBlock.stateFlag); tbl_regression.actionXbias = allMouse{i}.biasBlock.stateFlag.* tbl_regression.action;
    tbl_regression = tbl_regression(trialNum_bias{i,1},:); 
    tbl_regression(tbl_regression.context==3,:) = [];tbl_regression.context = [];
        
    tbl_ctxt.RT = tempRT; tbl_ctxt.stateFlag = allMouse{i}.biasBlock.stateFlag; 
    tbl_ctxt.onTime = tempAttributes(:,1); tbl_ctxt.onSpeed = tempAttributes(:,2); tbl_ctxt.onDist = tempAttributes(:,3);
    tbl_ctxt.totalTime = tempAttributes(:,4); tbl_ctxt.totalSpeed = tempAttributes(:,5); tbl_ctxt.totalDist = tempAttributes(:,6);

    [wheelDownSample, ~,selFlag] = loadWheelClusterData(mouseMega, i);
    templabel = label4(selFlag);  tbl_ctxt.label = templabel'; 
    
    tbl_ctxt = tbl_ctxt(trialNum_bias{i,1},:);
    tbl_ctxt(tbl_ctxt.context==3,:) = [];

    %temp = unique(tbl_regression.day);
    %tempDayFlag =  tbl_regression.day == temp(3); disp(sum(tempDayFlag))

    allTbl_bias{i} = tbl_ctxt; allTbl_regression{i} = tbl_regression;
    allWheel{i} = wheelDownSample(trialNum_bias{i,1},:);

    
    labelCell{i} = tbl_ctxt.label; labelCell{i}(labelCell{i}<=2) = 0; labelCell{i}(labelCell{i}>2) = 1;
    RTCell{i} = tbl_ctxt.RT; attributeCell{i} = table2array(tbl_ctxt(:,end-6:end-1));
end

%% PLOT 4.01 -- GET PERIOD OF BIAS VS UNBIASED STATE, for old method of bias block detection, transitioning out of bias
probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
outCell = mouseMega.objFun('binProbeByTrialFromLearningOnset',{[reinfThre probeThre],probeTrialBin});
learningOnset = cell2mat(outCell{1,3}); learningStart = learningOnset - 800; learningEnd = learningOnset + 800;
learningStart(learningStart<=100) = 100; 

nTrialBefAftBlock = 20;
trialNum_bias = cell(mouseMega.nMouse,3); allTbl_bias = {}; allTbl_regression = {}; allWheel = {};
RTcell = {};labelCell = {};attributeCell = {};
for i = 1:mouseMega.nMouse
    tempL = allMouse{i}.biasBlock.blockL.end; tempL_len = allMouse{i}.biasBlock.blockL.len; 
    tempR = allMouse{i}.biasBlock.blockR.end; tempR_len = allMouse{i}.biasBlock.blockR.len; 
    for j = 1:length(tempL)
        temp = min([tempL_len(j) nTrialBefAftBlock]);
        trialNum_bias{i,2} = cat(2,trialNum_bias{i,2},tempL(j)-temp:tempL(j)+temp);  
    end
    for j = 1:length(tempR)
        temp = min([tempR_len(j) nTrialBefAftBlock]);
        trialNum_bias{i,3} = cat(2,trialNum_bias{i,3},tempR(j)-temp:tempR(j)+temp);
    end
    trialNum_bias{i,2} = unique(trialNum_bias{i,2}); trialNum_bias{i,3} = unique(trialNum_bias{i,3});
    trialNum_bias{i,2}( trialNum_bias{i,2}<1) = []; trialNum_bias{i,3}( trialNum_bias{i,3}<1) = [];
    %trialNum_bias{i,2}(trialNum_bias{i,2}<learningStart(i) | trialNum_bias{i,2}>learningEnd(i)) = [];
    %trialNum_bias{i,3}(trialNum_bias{i,3}<learningStart(i) | trialNum_bias{i,3}>learningEnd(i)) = [];
    trialNum_bias{i,1} = cat(2, trialNum_bias{i,2}, trialNum_bias{i,3});


    [tbl_ctxt,tempRT,tempLabel,tempAttributes,~] = makeRegressorTable(mouseMega,i,label4,attributesOrig); 
    tbl_regression = tbl_ctxt;
    tbl_regression.bias = abs(allMouse{i}.biasBlock.stateFlag); tbl_regression.actionXbias = allMouse{i}.biasBlock.stateFlag.* tbl_regression.action;
    tbl_regression = tbl_regression(trialNum_bias{i,1},:); 
    tbl_regression(tbl_regression.context==3,:) = [];tbl_regression.context = [];
        
    tbl_ctxt.RT = tempRT; tbl_ctxt.stateFlag = allMouse{i}.biasBlock.stateFlag; 
    tbl_ctxt.onTime = tempAttributes(:,1); tbl_ctxt.onSpeed = tempAttributes(:,2); tbl_ctxt.onDist = tempAttributes(:,3);
    tbl_ctxt.totalTime = tempAttributes(:,4); tbl_ctxt.totalSpeed = tempAttributes(:,5); tbl_ctxt.totalDist = tempAttributes(:,6);

    [wheelDownSample, ~,selFlag] = loadWheelClusterData(mouseMega, i);
    templabel = label4(selFlag);  tbl_ctxt.label = templabel'; 
    
    tbl_ctxt = tbl_ctxt(trialNum_bias{i,1},:);
    tbl_ctxt(tbl_ctxt.context==3,:) = [];

    %temp = unique(tbl_regression.day);
    %tempDayFlag =  tbl_regression.day == temp(3); disp(sum(tempDayFlag))

    allTbl_bias{i} = tbl_ctxt; allTbl_regression{i} = tbl_regression;
    allWheel{i} = wheelDownSample(trialNum_bias{i,1},:);

    
    labelCell{i} = tbl_ctxt.label; labelCell{i}(labelCell{i}<=2) = 0; labelCell{i}(labelCell{i}>2) = 1;
    RTCell{i} = tbl_ctxt.RT; attributeCell{i} = table2array(tbl_ctxt(:,end-6:end-1));
end

%% PLOT 4.1 -- GET PERIOD OF BIAS VS UNBIASED STATE, for NEW method of bias block detection
probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
outCell = mouseMega.objFun('binProbeByTrialFromLearningOnset',{[reinfThre probeThre],probeTrialBin});
learningOnset = cell2mat(outCell{1,3}); learningStart = learningOnset - 800; learningEnd = learningOnset + 800;
learningStart(learningStart<=100) = 100; 
biasThre = 0.4; biasTrialThre = 80;

nTrialBefAftBlock = 20; nTrialAftBlock = 30; 
trialNum_bias = cell(mouseMega.nMouse,3); allTbl_bias = {}; allTbl_regression = {}; allWheel = {};
RTcell = {};labelCell = {};attributeCell = {};
trialDiff = [];
for i = 1:mouseMega.nMouse
    tempL = allMouse{i}.biasBlock.blockL.start; tempL_len = allMouse{i}.biasBlock.blockL.len; 
    tempR = allMouse{i}.biasBlock.blockR.start; tempR_len = allMouse{i}.biasBlock.blockR.len; 
    for j = 1:length(tempL)
        temp = min([tempL_len(j) nTrialBefAftBlock]);
        trialNum_bias{i,2} = cat(2,trialNum_bias{i,2},tempL(j)-nTrialBefAftBlock:tempL(j)-1); 
        tempIdx = tempL(j):allMouse{i}.biasBlock.blockL.end(j);
        tempBiasFlag = abs(allMouse{i}.behav.bias(tempIdx)) >= biasThre;tempIdx = tempIdx(tempBiasFlag);
        tempIdx(tempIdx-tempL(j) > biasTrialThre) = [];
        trialDiff(end+1) = mean(tempIdx) - tempL(j);
        trialNum_bias{i,2} = cat(2,trialNum_bias{i,2},tempIdx);
        %temp = min([tempL_len(j) nTrialBefAftBlock + nTrialAftBlock]);
        %trialNum_bias{i,2} = cat(2,trialNum_bias{i,2},tempL(j)-nTrialBefAftBlock:tempL(j)-1); 
        %trialNum_bias{i,2} = cat(2,trialNum_bias{i,2},tempL(j)+nTrialBefAftBlock:tempL(j)+temp); 
    end
    for j = 1:length(tempR)
        temp = min([tempR_len(j) nTrialBefAftBlock]);
        trialNum_bias{i,3} = cat(2,trialNum_bias{i,3},tempR(j)-nTrialBefAftBlock:tempR(j)-1);
        tempIdx = tempR(j):allMouse{i}.biasBlock.blockR.end(j);
        tempBiasFlag = abs(allMouse{i}.behav.bias(tempIdx)) >= biasThre;tempIdx = tempIdx(tempBiasFlag);
        tempIdx(tempIdx-tempR(j) > biasTrialThre) = [];
        trialDiff(end+1) = mean(tempIdx) - tempR(j);
        trialNum_bias{i,2} = cat(2,trialNum_bias{i,2},tempIdx);
        %temp = min([tempR_len(j) nTrialBefAftBlock + nTrialAftBlock]);
        %trialNum_bias{i,3} = cat(2,trialNum_bias{i,3},tempR(j)-nTrialBefAftBlock:tempR(j)-1);
        %trialNum_bias{i,3} = cat(2,trialNum_bias{i,3},tempR(j)+nTrialBefAftBlock:tempR(j)+temp);
    end
    trialNum_bias{i,2} = unique(trialNum_bias{i,2}); trialNum_bias{i,3} = unique(trialNum_bias{i,3});
    %trialNum_bias{i,2}( trialNum_bias{i,2}<1) = []; trialNum_bias{i,3}( trialNum_bias{i,3}<1) = [];
    trialNum_bias{i,2}(trialNum_bias{i,2}<learningStart(i) | trialNum_bias{i,2}>learningEnd(i)) = [];
    trialNum_bias{i,3}(trialNum_bias{i,3}<learningStart(i) | trialNum_bias{i,3}>learningEnd(i)) = [];
    trialNum_bias{i,1} = cat(2, trialNum_bias{i,2}, trialNum_bias{i,3});


    [tbl_ctxt,tempRT,tempLabel,tempAttributes,~] = makeRegressorTable(mouseMega,i,label4,attributesOrig); 
    tbl_regression = tbl_ctxt;
    %tbl_regression.bias = abs(allMouse{i}.biasBlock.stateFlag); tbl_regression.actionXbias = allMouse{i}.biasBlock.stateFlag.* tbl_regression.action; 
    tbl_regression.bias = abs(allMouse{i}.biasBlock.stateFlag)'; tbl_regression.actionXbias = allMouse{i}.biasBlock.stateFlag' .* tbl_regression.action; 
    tbl_regression = tbl_regression(trialNum_bias{i,1},:); 
    tbl_regression(tbl_regression.context==3,:) = [];tbl_regression.context = [];
        
    tbl_ctxt.RT = tempRT; %tbl_ctxt.stateFlag = allMouse{i}.biasBlock.stateFlag; 
    tbl_ctxt.stateFlag = allMouse{i}.biasBlock.stateFlag'; 
    tbl_ctxt.onTime = tempAttributes(:,1); tbl_ctxt.onSpeed = tempAttributes(:,2); tbl_ctxt.onDist = tempAttributes(:,3);
    tbl_ctxt.totalTime = tempAttributes(:,4); tbl_ctxt.totalSpeed = tempAttributes(:,5); tbl_ctxt.totalDist = tempAttributes(:,6);

    [wheelDownSample, ~,selFlag] = loadWheelClusterData(mouseMega, i);
    templabel = label4(selFlag);  tbl_ctxt.label = templabel'; 
    
    tbl_ctxt = tbl_ctxt(trialNum_bias{i,1},:);
    tbl_ctxt(tbl_ctxt.context==3,:) = [];

    %temp = unique(tbl_regression.day);
    %tempDayFlag =  tbl_regression.day == temp(3); disp(sum(tempDayFlag))

    allTbl_bias{i} = tbl_ctxt; allTbl_regression{i} = tbl_regression;
    allWheel{i} = wheelDownSample(trialNum_bias{i,1},:);

    
    labelCell{i} = tbl_ctxt.label; labelCell{i}(labelCell{i}<=2) = 0; labelCell{i}(labelCell{i}>2) = 1;
    RTCell{i} = tbl_ctxt.RT; attributeCell{i} = table2array(tbl_ctxt(:,end-6:end-1));
end
figure; histogram(trialDiff)
%%
plotBiasBlockDiff(allTbl_bias,'RT');
%%
plotBiasBlockDiff(allTbl_bias,'onTime');
plotBiasBlockDiff(allTbl_bias,'onSpeed');
plotBiasBlockDiff(allTbl_bias,'onDist');
%%
plotBiasBlockDiff(allTbl_bias,'totalTime');
plotBiasBlockDiff(allTbl_bias,'totalSpeed');
plotBiasBlockDiff(allTbl_bias,'totalDist');
%%
plotBiasBlockDiffDir(allTbl_bias,'RT',1);
plotBiasBlockDiffDir(allTbl_bias,'RT',-1);

%%
plotBiasBlockDiffDir(allTbl_bias,'onSpeed',1);
plotBiasBlockDiffDir(allTbl_bias,'onSpeed',-1);

%%

removeVarNames = {'rewardH3','rewardH2','rewardH1','day','trials','bias','choiceXbias','actionRate'};
tblLM_rsq = []; AIC = [];
for i = 1:(length(removeVarNames)+1)
    if i~=1
        for j = 1:length(allTbl_regression); allTbl_regression{j}.(removeVarNames{i-1}) = []; end
    end
    [tblLM_rsq_temp,~, ~, ~,AIC_temp, ~] = cellfun(@(a,b,c,d)(fitLM_allVar(a,b,c,d,false)), allTbl_regression...
        ,RTCell,labelCell,attributeCell, 'UniformOutput', false);
    temp = tblLM_rsq_temp;tblLM_rsq_temp = table();
    for k = 1:mouseMega.nMouse; tblLM_rsq_temp(:,k) = temp{k}; end; tblLM_rsq(:,:,i) = table2array(tblLM_rsq_temp);
    temp = AIC_temp;AIC_temp = table();
    for k = 1:mouseMega.nMouse; AIC_temp(:,k) = temp{k}; end; AIC(:,:,i) = table2array(AIC_temp);
end
%%
modelNames = {'+rewardH3','+rewardH2','rewardH1','+day','+trials','+bias','+choiceXbias','+actionRate','action-only'};
plotModelComparison(AIC,tblLM_rsq_temp.Properties.RowNames,modelNames,'AIC');
plotModelComparison(tblLM_rsq,tblLM_rsq_temp.Properties.RowNames,modelNames,'Var Exp');

selModel = 6; selModelFit = tblLM_rsq(:,:,selModel);
figure;  fn_plotBar(selModelFit); xticklabels(tblLM_rsq_temp.Properties.RowNames); ylabel('R-squared')

selModel = 4; selModelFit = tblLM_rsq(:,:,selModel);
figure;  fn_plotBar(selModelFit); xticklabels(tblLM_rsq_temp.Properties.RowNames); ylabel('R-squared')

%% PART 6.0 -- look at if reaction time difference (just between two actions) stay in the same direction
f1 = figure; f2 = figure;
for i = 1:6
    [~,tempRT,tempLabel,~,~] = makeRegressorTable(mouseMega,i,label4,attributesOrig); 

    aL = tempRT; aR = tempRT;
    b = allMouse{i}.behav.action;
    aL(b==2 | tempLabel==1 ) = nan; aR(b==1 | tempLabel==1) = nan;
    c1 = smoothdata(aL,'movmean',100);c2 = smoothdata(aR,'movmean',100);
    figure(f1);subplot(2,3,i); plot(c2-c1); hold on; plot(allMouse{i}.behav.bias); title(corr(c2-c1,allMouse{i}.behav.bias,'rows','complete'))
    figure(f2);subplot(2,3,i); hold on; plot([1,length(c1)],[0 0],'Color',[0.8 0.8 0.8]); plot(c2-c1,'Color',matlabColors(1)); 
end



%%
biasTransL = {}; biasTransR = {};
biasTransRTL = {}; biasTransRTR = {};
for i = 1:length(allMouse)
    biasTransL{i} = {};  biasTransR{i} = {}; biasTransRTL{i} = [];  biasTransRTR{i} = [];
    [wheelDownSample_all, RT_all,selFlag] = loadWheelClusterData(mouseMega, i);
    tempLabel = label4(selFlag);

    for j = 1:length(allMouse{i}.biasBlock.transID)



        if strcmp(allMouse{i}.biasBlock.transID{j},'U2L') || strcmp(allMouse{i}.biasBlock.transID{j},'U2R') ||...
            strcmp(allMouse{i}.biasBlock.transID{j},'L2R') || strcmp(allMouse{i}.biasBlock.transID{j},'R2L')
            tempTrans = allMouse{i}.biasBlock.trans(j,:);
            if tempTrans(2)-tempTrans(1) <= 30
                tempUnbiasIdx = tempTrans(1):tempTrans(2);
            else 
                tempUnbiasIdx = tempTrans(2)-30:tempTrans(2)-1;
            end

         
            if tempTrans(4) - tempTrans(3) <= 70
                tempLen = floor((tempTrans(4) - tempTrans(3)+1)/3);
                tempBiasIdx = tempTrans(3)+tempLen:tempTrans(3)+tempLen*2;
            else
                tempBiasIdx = tempTrans(3)+21:tempTrans(3)+50;

            end
            RT_unbias = nanmean(RT_all(fn_idx2logical(tempUnbiasIdx,length(tempLabel))' & (tempLabel<=2)));
            RT_bias = nanmean(RT_all(fn_idx2logical(tempBiasIdx,length(tempLabel))' & (tempLabel<=2)));



            if strcmp(allMouse{i}.biasBlock.transID{j},'U2L') || strcmp(allMouse{i}.biasBlock.transID{j},'R2L')
               
                
                
                biasTransL{i}(end+1,:) = {tempUnbiasIdx, tempBiasIdx};
                biasTransRTL{i}(end+1,:) = [RT_unbias,RT_bias];
            elseif strcmp(allMouse{i}.biasBlock.transID{j},'U2R') || strcmp(allMouse{i}.biasBlock.transID{j},'L2R')
                biasTransR{i}(end+1,:) = {tempUnbiasIdx, tempBiasIdx};
                biasTransRTR{i}(end+1,:) = [RT_unbias,RT_bias];

            end

            


        end

    end

end
a = cellfun(@(x)(mean(x,1)),biasTransRTL,'UniformOutput',false);
b = cellfun(@(x)(mean(x,1)),biasTransRTR,'UniformOutput',false);
%% FUNCTIONS


function plotBiasBlockDiff(allTbl_ctxt,varName)
reinfLabelCount = zeros(length(allTbl_ctxt),2); probeLabelCount = zeros(length(allTbl_ctxt),2);
reinfRT = zeros(length(allTbl_ctxt),2); probeRT = zeros(length(allTbl_ctxt),2);
for i = 1:length(allTbl_ctxt)
    tempProbeL = allTbl_ctxt{i}.label(allTbl_ctxt{i}.stateFlag ~=0 & allTbl_ctxt{i}.action ==1);
    tempProbeR = allTbl_ctxt{i}.label(allTbl_ctxt{i}.stateFlag ~=0 & allTbl_ctxt{i}.action ==2); 

    probeLabelCount(i,:) = [nanmean(tempProbeL<=2) nanmean(tempProbeR<=2)];


    tempReinfL = allTbl_ctxt{i}.label(allTbl_ctxt{i}.stateFlag ==0 & allTbl_ctxt{i}.action ==1); 
    tempReinfR = allTbl_ctxt{i}.label(allTbl_ctxt{i}.stateFlag ==0 & allTbl_ctxt{i}.action ==2); 

    reinfLabelCount(i,:) = [nanmean(tempReinfL<=2) nanmean(tempReinfR<=2)];

    probeRT(i,1) = nanmean(allTbl_ctxt{i}.(varName)(allTbl_ctxt{i}.stateFlag ~=0 & allTbl_ctxt{i}.action ==1));
    probeRT(i,2) = nanmean(allTbl_ctxt{i}.(varName)(allTbl_ctxt{i}.stateFlag ~=0 & allTbl_ctxt{i}.action ==2));
    reinfRT(i,1) = nanmean(allTbl_ctxt{i}.(varName)(allTbl_ctxt{i}.stateFlag ==0 & allTbl_ctxt{i}.action ==1));
    reinfRT(i,2) = nanmean(allTbl_ctxt{i}.(varName)(allTbl_ctxt{i}.stateFlag ==0 & allTbl_ctxt{i}.action ==2));
end

figure; subplot(4,3,1); fn_plotBarPaired({nanmean(reinfLabelCount,2),nanmean(probeLabelCount,2)}); xticklabels({'unbiased','biased'}); ylabel('Proportion of fast cluster');
subplot(4,3,2); fn_plotBarPaired({reinfLabelCount(:,1),probeLabelCount(:,1)}); xticklabels({'unbiased','biased'}); ylabel('Proportion of fast cluster in Left');
subplot(4,3,3); fn_plotBarPaired({reinfLabelCount(:,2),probeLabelCount(:,2)}); xticklabels({'unbiased','biased'}); ylabel('Proportion of fast cluster in Right');
subplot(4,3,7); fn_plotBarPaired({nanmean(reinfRT,2),nanmean(probeRT,2)}); xticklabels({'unbiased','biased'}); ylabel([varName ' of all clusters']);
subplot(4,3,8); fn_plotBarPaired({reinfRT(:,1),probeRT(:,1)}); xticklabels({'unbiased','biased'}); ylabel([varName ' of Left Action']);
subplot(4,3,9); fn_plotBarPaired({reinfRT(:,2),probeRT(:,2)}); xticklabels({'unbiased','biased'}); ylabel([varName ' of Right Action']);
% LOOK AT FAST 1 VS. FAST 2 CLUSTERS
reinfLabelCount = zeros(length(allTbl_ctxt),3); probeLabelCount = zeros(length(allTbl_ctxt),3);
reinfRT = zeros(length(allTbl_ctxt),3); probeRT = zeros(length(allTbl_ctxt),3);

preferredDirection = [2 1 1 1 1 1];

for i = 1:length(allTbl_ctxt)
    probeRT(i,2) = (nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ~=0 | allTbl_ctxt{i}.label ==1) & allTbl_ctxt{i}.action ==1)) + ...
        nanmean( allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ~=0 | allTbl_ctxt{i}.label ==1) & allTbl_ctxt{i}.action ==2 )))/2;
    probeRT(i,3) = (nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ~=0 | allTbl_ctxt{i}.label ==2) & allTbl_ctxt{i}.action ==1)) + ...
        nanmean( allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ~=0 | allTbl_ctxt{i}.label ==2) & allTbl_ctxt{i}.action ==2 )))/2;
    probeRT(i,1) = (nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ~=0 | allTbl_ctxt{i}.label <=2) & allTbl_ctxt{i}.action ==1)) + ...
        nanmean( allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ~=0 | allTbl_ctxt{i}.label <=2) & allTbl_ctxt{i}.action ==2 )))/2;
    
    probeRT(i,4) = nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ~=0 | allTbl_ctxt{i}.label <=2) & allTbl_ctxt{i}.action == preferredDirection(i)));
    probeRT(i,5) = nanmean( allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ~=0 | allTbl_ctxt{i}.label <=2) & allTbl_ctxt{i}.action == (3-preferredDirection(i)) ));

    reinfRT(i,2) = (nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==0 | allTbl_ctxt{i}.label ==1) & allTbl_ctxt{i}.action ==1)) + ...
        nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==0 | allTbl_ctxt{i}.label ==1) & allTbl_ctxt{i}.action ==2)))/2;
    reinfRT(i,3) = (nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==0 | allTbl_ctxt{i}.label ==2) & allTbl_ctxt{i}.action ==1)) + ...
        nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==0 | allTbl_ctxt{i}.label ==2) & allTbl_ctxt{i}.action ==2)))/2;
    reinfRT(i,1) = (nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==0 | allTbl_ctxt{i}.label <=2) & allTbl_ctxt{i}.action ==1)) + ...
        nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==0 | allTbl_ctxt{i}.label <=2) & allTbl_ctxt{i}.action ==2)))/2;

    reinfRT(i,4) = nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==0 | allTbl_ctxt{i}.label <=2) & allTbl_ctxt{i}.action == preferredDirection(i) ));
    reinfRT(i,5) = nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==0 | allTbl_ctxt{i}.label <=2) & allTbl_ctxt{i}.action == (3-preferredDirection(i)) ));

    tempProbeL = allTbl_ctxt{i}.label(allTbl_ctxt{i}.stateFlag ~=0 & allTbl_ctxt{i}.action ==1);
    tempProbeR = allTbl_ctxt{i}.label(allTbl_ctxt{i}.stateFlag ~=0 & allTbl_ctxt{i}.action ==2); 
    tempL1 = nanmean(tempProbeL==1); tempR1 = nanmean(tempProbeR==1); 
    tempL2 = nanmean(tempProbeL==2); tempR2 = nanmean(tempProbeR==2); tempL = tempL1/(tempL1+tempL2); tempR = tempR1/(tempR1+tempR2);
    probeLabelCount(i,:) = [(tempL+tempR)/2 tempL tempR];

    tempReinfL = allTbl_ctxt{i}.label(allTbl_ctxt{i}.stateFlag ==0 & allTbl_ctxt{i}.action ==1); 
    tempReinfR = allTbl_ctxt{i}.label(allTbl_ctxt{i}.stateFlag ==0 & allTbl_ctxt{i}.action ==2); 
    tempL1 = nansum(tempReinfL==1); tempR1 = nansum(tempReinfR==1); 
    tempL2 = nansum(tempReinfL==2); tempR2 = nansum(tempReinfR==2); tempL = tempL1/(tempL1+tempL2); tempR = tempR1/(tempR1+tempR2);
    reinfLabelCount(i,:) = [(tempL+tempR)/2 tempL tempR];

    %tempProbe = allTbl_ctxt{i}.label(allTbl_ctxt{i}.context ==2); probeLabelCount(i,:) = [sum(tempProbe==1) sum(tempProbe==2)];
    %tempReinf = allTbl_ctxt{i}.label(allTbl_ctxt{i}.context ==1); reinfLabelCount(i,:) = [sum(tempReinf==1) sum(tempReinf==2)];
end

subplot(4,3,4); fn_plotBarPaired({reinfLabelCount(:,1),probeLabelCount(:,1)}); xticklabels({'unbiased','biased'}); ylabel('Proportion of type 1 fast cluster');

subplot(4,3,5); fn_plotBarPaired({reinfLabelCount(:,2),probeLabelCount(:,2)}); xticklabels({'unbiased','biased'}); ylabel('Proportion of type 1 fast cluster, Left');
subplot(4,3,6); fn_plotBarPaired({reinfLabelCount(:,3),probeLabelCount(:,3)}); xticklabels({'unbiased','biased'}); ylabel('Proportion of type 1 fast cluster, Right');

subplot(4,3,10); fn_plotBarPaired({reinfRT(:,1),probeRT(:,1)}); xticklabels({'unbiased','biased'}); ylabel([varName ' of fast clusters']);
subplot(4,3,11); fn_plotBarPaired({reinfRT(:,4),probeRT(:,4)}); xticklabels({'unbiased','biased'}); ylabel([varName ' of fast clusters, Left']);
subplot(4,3,12); fn_plotBarPaired({reinfRT(:,5),probeRT(:,5)}); xticklabels({'unbiased','biased'}); ylabel([varName ' of fast clusters, Right']);


figure;
subplot(4,2,1); fn_plotBarPaired({reinfRT(:,2),probeRT(:,2)}); xticklabels({'unbiased','biased'}); ylabel([varName ' of fast cluster 1']);
subplot(4,2,2); fn_plotBarPaired({reinfRT(:,3),probeRT(:,3)}); xticklabels({'unbiased','biased'}); ylabel([varName ' of fast cluster 2']);

subplot(4,2,3); fn_plotBarPaired(probeRT(:,2:3)'); xticklabels({'cluster 1','cluster 2'}); ylabel([varName ' of biased']);
subplot(4,2,4); fn_plotBarPaired(reinfRT(:,2:3)'); xticklabels({'cluster 1','cluster 2'}); ylabel([varName ' of unbiased']);

subplot(4,2,5); fn_plotBarPaired(probeRT(:,4:5)'); xticklabels({'preferred','unpreferred'}); ylabel([varName ' of biased']);
subplot(4,2,6); fn_plotBarPaired(reinfRT(:,4:5)'); xticklabels({'preferred','unpreferred'}); ylabel([varName ' of unbiased']);

for k = 1:length(preferredDirection)
    if preferredDirection(k)==2
        reinfLabelCount(:,2:3) = reinfLabelCount(:,[3 2]); 
        probeLabelCount(:,2:3) = probeLabelCount(:,[3 2]); 
    end
end
subplot(4,2,7); fn_plotBarPaired(probeLabelCount(:,2:3)'); xticklabels({'preferred','unpreferred'}); ylabel('Proportion of type 1 fast cluster, biased');
subplot(4,2,8); fn_plotBarPaired(reinfLabelCount(:,2:3)'); xticklabels({'preferred','unpreferred'}); ylabel('Proportion of type 1 fast cluster, unbiased');

end

function plotBiasBlockDiffDir(allTbl_ctxt,varName,direction)
reinfLabelCount = zeros(length(allTbl_ctxt),2); probeLabelCount = zeros(length(allTbl_ctxt),2);
reinfRT = zeros(length(allTbl_ctxt),2); probeRT = zeros(length(allTbl_ctxt),2);
for i = 1:length(allTbl_ctxt)
    tempProbeL = allTbl_ctxt{i}.label(allTbl_ctxt{i}.stateFlag ==direction & allTbl_ctxt{i}.action ==1);
    tempProbeR = allTbl_ctxt{i}.label(allTbl_ctxt{i}.stateFlag ==direction & allTbl_ctxt{i}.action ==2); 

    probeLabelCount(i,:) = [nanmean(tempProbeL<=2) nanmean(tempProbeR<=2)];


    tempReinfL = allTbl_ctxt{i}.label(allTbl_ctxt{i}.stateFlag ==0 & allTbl_ctxt{i}.action ==1); 
    tempReinfR = allTbl_ctxt{i}.label(allTbl_ctxt{i}.stateFlag ==0 & allTbl_ctxt{i}.action ==2); 

    reinfLabelCount(i,:) = [nanmean(tempReinfL<=2) nanmean(tempReinfR<=2)];

    probeRT(i,1) = nanmean(allTbl_ctxt{i}.(varName)(allTbl_ctxt{i}.stateFlag ==direction & allTbl_ctxt{i}.action ==1));
    probeRT(i,2) = nanmean(allTbl_ctxt{i}.(varName)(allTbl_ctxt{i}.stateFlag ==direction & allTbl_ctxt{i}.action ==2));
    reinfRT(i,1) = nanmean(allTbl_ctxt{i}.(varName)(allTbl_ctxt{i}.stateFlag ==0 & allTbl_ctxt{i}.action ==1));
    reinfRT(i,2) = nanmean(allTbl_ctxt{i}.(varName)(allTbl_ctxt{i}.stateFlag ==0 & allTbl_ctxt{i}.action ==2));
end

figure; subplot(4,3,1); fn_plotBarPaired({nanmean(reinfLabelCount,2),nanmean(probeLabelCount,2)}); xticklabels({'unbiased','biased'}); ylabel('Proportion of fast cluster');
subplot(4,3,2); fn_plotBarPaired({reinfLabelCount(:,1),probeLabelCount(:,1)}); xticklabels({'unbiased','biased'}); ylabel('Proportion of fast cluster in Left');
subplot(4,3,3); fn_plotBarPaired({reinfLabelCount(:,2),probeLabelCount(:,2)}); xticklabels({'unbiased','biased'}); ylabel('Proportion of fast cluster in Right');
subplot(4,3,7); fn_plotBarPaired({nanmean(reinfRT,2),nanmean(probeRT,2)}); xticklabels({'unbiased','biased'}); ylabel([varName ' of all clusters']);
subplot(4,3,8); fn_plotBarPaired({reinfRT(:,1),probeRT(:,1)}); xticklabels({'unbiased','biased'}); ylabel([varName ' of Left Action']);
subplot(4,3,9); fn_plotBarPaired({reinfRT(:,2),probeRT(:,2)}); xticklabels({'unbiased','biased'}); ylabel([varName ' of Right Action']);
% LOOK AT FAST 1 VS. FAST 2 CLUSTERS
reinfLabelCount = zeros(length(allTbl_ctxt),3); probeLabelCount = zeros(length(allTbl_ctxt),3);
reinfRT = zeros(length(allTbl_ctxt),3); probeRT = zeros(length(allTbl_ctxt),3);
for i = 1:length(allTbl_ctxt)
    probeRT(i,2) = (nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==direction | allTbl_ctxt{i}.label ==1) & allTbl_ctxt{i}.action ==1)) + ...
        nanmean( allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==direction | allTbl_ctxt{i}.label ==1) & allTbl_ctxt{i}.action ==2 )))/2;
    probeRT(i,3) = (nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==direction | allTbl_ctxt{i}.label ==2) & allTbl_ctxt{i}.action ==1)) + ...
        nanmean( allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==direction | allTbl_ctxt{i}.label ==2) & allTbl_ctxt{i}.action ==2 )))/2;
    probeRT(i,1) = (nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==direction | allTbl_ctxt{i}.label <=2) & allTbl_ctxt{i}.action ==1)) + ...
        nanmean( allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==direction | allTbl_ctxt{i}.label <=2) & allTbl_ctxt{i}.action ==2 )))/2;
    
    probeRT(i,4) = nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==direction | allTbl_ctxt{i}.label <=2) & allTbl_ctxt{i}.action ==1));
    probeRT(i,5) = nanmean( allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==direction | allTbl_ctxt{i}.label <=2) & allTbl_ctxt{i}.action ==2 ));

    reinfRT(i,2) = (nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==0 | allTbl_ctxt{i}.label ==1) & allTbl_ctxt{i}.action ==1)) + ...
        nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==0 | allTbl_ctxt{i}.label ==1) & allTbl_ctxt{i}.action ==2)))/2;
    reinfRT(i,3) = (nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==0 | allTbl_ctxt{i}.label ==2) & allTbl_ctxt{i}.action ==1)) + ...
        nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==0 | allTbl_ctxt{i}.label ==2) & allTbl_ctxt{i}.action ==2)))/2;
    reinfRT(i,1) = (nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==0 | allTbl_ctxt{i}.label <=2) & allTbl_ctxt{i}.action ==1)) + ...
        nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==0 | allTbl_ctxt{i}.label <=2) & allTbl_ctxt{i}.action ==2)))/2;

    reinfRT(i,4) = nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==0 | allTbl_ctxt{i}.label <=2) & allTbl_ctxt{i}.action ==1));
    reinfRT(i,5) = nanmean(allTbl_ctxt{i}.(varName)( (allTbl_ctxt{i}.stateFlag ==0 | allTbl_ctxt{i}.label <=2) & allTbl_ctxt{i}.action ==2));

    tempProbeL = allTbl_ctxt{i}.label(allTbl_ctxt{i}.stateFlag ==direction & allTbl_ctxt{i}.action ==1);
    tempProbeR = allTbl_ctxt{i}.label(allTbl_ctxt{i}.stateFlag ==direction & allTbl_ctxt{i}.action ==2); 
    tempL1 = nanmean(tempProbeL==1); tempR1 = nanmean(tempProbeR==1); 
    tempL2 = nanmean(tempProbeL==2); tempR2 = nanmean(tempProbeR==2); tempL = tempL1/(tempL1+tempL2); tempR = tempR1/(tempR1+tempR2);
    probeLabelCount(i,:) = [(tempL+tempR)/2 tempL tempR];

    tempReinfL = allTbl_ctxt{i}.label(allTbl_ctxt{i}.stateFlag ==0 & allTbl_ctxt{i}.action ==1); 
    tempReinfR = allTbl_ctxt{i}.label(allTbl_ctxt{i}.stateFlag ==0 & allTbl_ctxt{i}.action ==2); 
    tempL1 = nansum(tempReinfL==1); tempR1 = nansum(tempReinfR==1); 
    tempL2 = nansum(tempReinfL==2); tempR2 = nansum(tempReinfR==2); tempL = tempL1/(tempL1+tempL2); tempR = tempR1/(tempR1+tempR2);
    reinfLabelCount(i,:) = [(tempL+tempR)/2 tempL tempR];

    %tempProbe = allTbl_ctxt{i}.label(allTbl_ctxt{i}.context ==2); probeLabelCount(i,:) = [sum(tempProbe==1) sum(tempProbe==2)];
    %tempReinf = allTbl_ctxt{i}.label(allTbl_ctxt{i}.context ==1); reinfLabelCount(i,:) = [sum(tempReinf==1) sum(tempReinf==2)];
end

subplot(4,3,4); fn_plotBarPaired({reinfLabelCount(:,1),probeLabelCount(:,1)}); xticklabels({'unbiased','biased'}); ylabel('Proportion of type 1 fast cluster');

subplot(4,3,5); fn_plotBarPaired({reinfLabelCount(:,2),probeLabelCount(:,2)}); xticklabels({'unbiased','biased'}); ylabel('Proportion of type 1 fast cluster, Left');
subplot(4,3,6); fn_plotBarPaired({reinfLabelCount(:,3),probeLabelCount(:,3)}); xticklabels({'unbiased','biased'}); ylabel('Proportion of type 1 fast cluster, Right');

subplot(4,3,10); fn_plotBarPaired({reinfRT(:,1),probeRT(:,1)}); xticklabels({'unbiased','biased'}); ylabel([varName ' of fast clusters']);
subplot(4,3,11); fn_plotBarPaired({reinfRT(:,4),probeRT(:,4)}); xticklabels({'unbiased','biased'}); ylabel([varName ' of fast clusters, Left']);
subplot(4,3,12); fn_plotBarPaired({reinfRT(:,5),probeRT(:,5)}); xticklabels({'unbiased','biased'}); ylabel([varName ' of fast clusters, Right']);


figure;
subplot(4,2,1); fn_plotBarPaired({reinfRT(:,2),probeRT(:,2)}); xticklabels({'unbiased','biased'}); ylabel([varName ' of fast cluster 1']);
subplot(4,2,2); fn_plotBarPaired({reinfRT(:,3),probeRT(:,3)}); xticklabels({'unbiased','biased'}); ylabel([varName ' of fast cluster 2']);

subplot(4,2,3); fn_plotBarPaired(probeRT(:,2:3)'); xticklabels({'cluster 1','cluster 2'}); ylabel([varName ' of biased']);
subplot(4,2,4); fn_plotBarPaired(reinfRT(:,2:3)'); xticklabels({'cluster 1','cluster 2'}); ylabel([varName ' of unbiased']);

subplot(4,2,5); fn_plotBarPaired(probeRT(:,4:5)'); xticklabels({'left','right'}); ylabel([varName ' of biased']);
subplot(4,2,6); fn_plotBarPaired(reinfRT(:,4:5)'); xticklabels({'left','right'}); ylabel([varName ' of unbiased']);

subplot(4,2,7); fn_plotBarPaired(probeLabelCount(:,2:3)'); xticklabels({'left','right'}); ylabel('Proportion of type 1 fast cluster, biased');
subplot(4,2,8); fn_plotBarPaired(reinfLabelCount(:,2:3)'); xticklabels({'left','right'}); ylabel('Proportion of type 1 fast cluster, unbiased');

end

function plotModelComparison(metric,targVarNames,modelNames,metricName)
% Look at the improvement for each animal
%{
for i = 1:size(AIC,1)
    figure;
    for j = 1:size(AIC,2)
        subplot(2,3,j)
        plot(squeeze(AIC(i,j,:)))
    end
end
%}
metric = metric - repmat(metric(:,:,end),[1 1 size(metric,3)]);

figure; 
for i = 1:size(metric,1)
    subplot(2,5,i); hold on; 
    fn_plotMeanErrorbar(1:size(metric,3),squeeze(metric(i,:,:)),matlabColors(1,0.8),matlabColors(1,0.2),{},{});
    xticks(1:9);xticklabels(modelNames); ylabel(['Delta ' metricName]);
    title(targVarNames{i})
end

end

function [tblLM_rsq,tblLM_coeff, tblLM_coeffP, tblLM_uniContr,BIC, pred] = fitLM_allVar(tbl,RT, label, attributesScaled,uniqueFlag)
    tbl_fast = tbl(label==0,:); RT_fast = RT(label==0); 
    tbl_slow = tbl(label==1,:); RT_slow = RT(label==1); 
    
    rowNames = {'RT','RTfast','RTslow','label','onTime','onSpeed','onDist','totalTime','totalSpeed','totalDist'};

    tblCell = {tbl,tbl_fast,tbl_slow,tbl,tbl,tbl,tbl,tbl,tbl,tbl};
    responseCell = {RT,RT_fast,RT_slow,label,attributesScaled(:,1),attributesScaled(:,2),attributesScaled(:,3)...
        ,attributesScaled(:,4),attributesScaled(:,5),attributesScaled(:,6)};

    [tblLM_rsq,tblLM_coeff, tblLM_coeffP, tblLM_uniContr,BIC, pred] = fitLM(tblCell,responseCell,rowNames,uniqueFlag);
end

function [tblLM_rsq,tblLM_coeff, tblLM_coeffP, tblLM_uniContr,BIC, pred] = fitLM(tbl,targetVarCell,targetVarCellName,uniqueFlag)
    tblLM_rsq = table(); tblLM_coeff = table(); tblLM_coeffP = table(); pred = {};
    tblLM_uniContr = table(); BIC = table();

    for i = 1:length(targetVarCell)
        tempTbl = tbl{i}; tempTbl.(targetVarCellName{i}) = targetVarCell{i};
        [tblLM_rsq{i,1}, tblLM_coeff{i,:}, tblLM_coeffP{i,:}, BIC{i,1}, ~, pred{i}] = getLM(tempTbl, []);
        if uniqueFlag
            uniContr = fitLM_uniqueContribution(tempTbl,tblLM_rsq{i,1});
            tblLM_uniContr{i,:} = uniContr;
        end
    end

    regressorNames = ['Intercept',tbl{1}.Properties.VariableNames];
    tblLM_coeff.Properties.VariableNames = regressorNames;
    tblLM_coeffP.Properties.VariableNames = regressorNames;
    
    tblLM_rsq.Properties.RowNames = targetVarCellName;
    tblLM_coeff.Properties.RowNames = targetVarCellName;
    tblLM_coeffP.Properties.RowNames = targetVarCellName;
    if uniqueFlag
        tblLM_uniContr.Properties.VariableNames = tbl{1}.Properties.VariableNames;
        tblLM_uniContr.Properties.RowNames = targetVarCellName;
    end
end

function [uniContr] = fitLM_uniqueContribution(tempTbl,rsq_full)
    uniContr = [];
    allRegNames = tempTbl.Properties.VariableNames;
    % Fit individual model with each variable removed
    for i = 1:length(allRegNames)-1
        [rsq, ~, ~, ~,~,~] = getLM(tempTbl, allRegNames{i}); 
        uniContr(i) = rsq_full - rsq; 
    end
end



function [rsq, coeffEst, coeffP, BIC, lm, pred] = getLM(tbl,removeName)
    tbl(:,1:end-1) = normalize(tbl(:,1:end-1));
    if isempty(removeName)
        lm = fitlm(tbl,'ResponseVar',tbl.Properties.VariableNames{end});
        rsq = lm.Rsquared.Adjusted; coeffEst = lm.Coefficients.Estimate'; 
        coeffP = lm.Coefficients.pValue'; pred = lm.Fitted; BIC = lm.ModelCriterion.AIC;
    else
        [rsq,coeffEst,coeffP,lm] = fn_fitlmRemoveRegressor(tbl,removeName);
        pred = lm.Fitted;BIC = lm.ModelCriterion.AIC;
    end
end

function [rsq, coeffEst, coeffP,lm] = getGLM(tbl,removeName)
    tbl(:,1:end-1) = normalize(tbl(:,1:end-1));
    
   if isempty(removeName)
        lm = fitglm(tbl,'Link','logit','Distribution','Binomial','ResponseVar',tbl.Properties.VariableNames{end});
        rsq = lm.Rsquared.Ordinary; coeffEst = lm.Coefficients.Estimate'; 
        coeffP = lm.Coefficients.pValue';
    else
        [rsq,coeffEst,coeffP,lm] = fn_fitglmRemoveRegressor(tbl,removeName);
   end
end

function plotBiasRT(RT,bias,label)
idx_biased = bias >= 0.4 | bias < -0.4; idx_unbiased = bias < 0.2 & bias >= -0.2;
label_biased = label(idx_biased); label_unbiased = label(idx_unbiased);
k = 4; labelCount_biased = zeros(1,k); labelCount_unbiased = zeros(1,k);
for i = 1:k
    labelCount_biased(i) = sum(label_biased==i);
    labelCount_unbiased(i) = sum(label_unbiased==i);
end
labelCount_biased = labelCount_biased / sum(labelCount_biased); labelCount_unbiased = labelCount_unbiased / sum(labelCount_unbiased);

figure; subplot(2,2,1); bar([1 2],cat(1,labelCount_biased,labelCount_unbiased),'stacked'); title('Proportion of all clusters')
xticks([1 2]); xticklabels({'biased','unbiased'})

% only plot the bars of the 2nd 
labelCount_biased = labelCount_biased(1:2) / sum(labelCount_biased(1:2)); labelCount_unbiased = labelCount_unbiased(1:2) / sum(labelCount_unbiased(1:2));
subplot(2,2,2); bar([1 2],cat(1,labelCount_biased,labelCount_unbiased),'stacked'); title('Proportion of first two clusters')
xticks([1 2]); xticklabels({'biased','unbiased'})

RT_biased = RT (label' <=2 & idx_biased); RT_unbiased = RT (label' <=2 & idx_unbiased);
subplot(2,2,3); hold on; cdfplot(RT_biased); cdfplot(RT_unbiased); legend({'biased','unbiased'}); xlabel('RT'); title('Reaction time of cluster 1 and 2 trials')

RT_biased = RT (label' >2 & idx_biased); RT_unbiased = RT (label' >2 & idx_unbiased);
subplot(2,2,4); hold on; cdfplot(RT_biased); cdfplot(RT_unbiased); legend({'biased','unbiased'}); xlabel('RT'); title('Reaction time of cluster 1 and 2 trials')
end

function plotActionRT(RT,action,label,barLabel)
idxL = (action == 1); idxR = (action ==2);
labelL = label(idxL); labelR = label(idxR);
k = 4; labelCountL = zeros(1,k); labelCountR = zeros(1,k);
for i = 1:k
    labelCountL(i) = sum(labelL==i);
    labelCountR(i) = sum(labelR==i);
end
labelCountL = labelCountL / sum(labelCountL); labelCountR = labelCountR / sum(labelCountR);

figure; subplot(2,2,1); bar([1 2],cat(1,labelCountL,labelCountR),'stacked'); title('Proportion of all clusters')
xticks([1 2]); xticklabels(barLabel)

% only plot the bars of the 2nd 
labelCountL = labelCountL(1:2) / sum(labelCountL(1:2)); labelCountR = labelCountR(1:2) / sum(labelCountR(1:2));
subplot(2,2,2); bar([1 2],cat(1,labelCountL,labelCountR),'stacked'); title('Proportion of first two clusters')
xticks([1 2]); xticklabels(barLabel)

RTL = RT (label' <=2 & idxL); RTR = RT (label' <=2 & idxR);
subplot(2,2,3); hold on; cdfplot(RTL); cdfplot(RTR); legend(barLabel); xlabel('RT'); title('Reaction time of cluster 1 and 2 trials')

RTL = RT (label' >2 & idxL); RTR = RT (label' >2 & idxR);
subplot(2,2,4); hold on; cdfplot(RTL); cdfplot(RTR); legend(barLabel); xlabel('RT'); title('Reaction time of cluster 3 and 4 trials')
end


function [probeData,maxAlignPoint] = attachNan(probeData, alignPoint)
    attachDim = 1; 
    maxAlignPoint = max(alignPoint);
    tempFieldNames = fieldnames(probeData);
    for i = 1:length(tempFieldNames)
        tempField = probeData.(tempFieldNames{i});
        if isnumeric(tempField)
            tempMatSize = size(tempField);
            tempMatSize(attachDim) = tempMatSize(attachDim) + maxAlignPoint;
            tempMat = nan(tempMatSize);
            for j = 1:length(alignPoint)
                startPoint = maxAlignPoint-alignPoint(j)+1;
                tempMat(startPoint:startPoint+size(tempField,attachDim)-1,j) = tempField(:,j);

            end
            probeData.(tempFieldNames{i}) = tempMat;
        end        
    end
end

function saveWheelTrajAttr(allMouse, mice)
[wheelDownSample,badFlag] = getDownSample(allMouse,[]);
downSampleRate = 50;
for mouse = 1:length(allMouse)
    tempWheelDownSample = wheelDownSample{mouse};
    nTrials = size(tempWheelDownSample,1);
    
    tempWheelChange = [zeros(size(tempWheelDownSample,1),1) diff(tempWheelDownSample,1,2)];
    wheelOneBlockIdx = false(nTrials,1); wheelMultiBlockIdx = false(nTrials,1); tempOnsetIdx = zeros(nTrials,1);
    
    onsetDist = []; onsetSpeed = [];totalDist = []; onsetTime = []; totalTime = []; totalSpeed = [];  onsetFrame = [];
    for i = 1:nTrials
        % find one-block or multiple block movements
        [onsetIdx, offIdx,blockIdx] = fn_getBlockOnOff(tempWheelChange(i,:)~=0);
        % correct for short blocks
        tempSelFlag = [];
        for j = 1:length(blockIdx); tempSelFlag(j) = sum(abs(tempWheelChange(i,blockIdx{j}))); end
        moveBlockThre = 0.2; onsetIdx(tempSelFlag< moveBlockThre) = [];
        offIdx(tempSelFlag< moveBlockThre) = []; blockIdx(tempSelFlag< moveBlockThre) = [];
        if ~isempty(blockIdx)
            onsetDist(i) = sum(abs(tempWheelChange(i,blockIdx{1}))); 
            onsetSpeed(i) = onsetDist(i) / length(blockIdx{1}) * 1000/downSampleRate; % convert to ms for speed
            totalDist(i) = sum(abs(tempWheelChange(i,:))) ; 
            totalTime(i) = (offIdx(end) - onsetIdx(1)) * downSampleRate; % convert to ms
            totalSpeed(i) = totalDist(i)/totalTime(i) * 1000;
            onsetTime(i) = onsetIdx(1) * downSampleRate; onsetFrame(i) = onsetIdx(1);
        end 
            
        
        if length(onsetIdx)>1; wheelMultiBlockIdx(i) = 1;
        else; wheelOneBlockIdx(i) = 1; end      
        % align and do clustering
        %tempOnsetIdx(i) = find(abs(tempWheelDownSample(i,:))~=0,1); 
    end
    %[tempWheelDownSample_aligned,maxAlignPoint] = fn_align2idx(tempWheelDownSample, tempOnsetIdx,'fill','startEndValue');
    [tempWheelDownSample_aligned,maxAlignPoint] = fn_align2idx(tempWheelDownSample, onsetFrame,'fill','startEndValue');
    
    figure; subplot(2,3,1); histogram(onsetTime);title('Onset Time (ms)');xlim([0 1000]) %1. onset time
    subplot(2,3,2); histogram(onsetSpeed); title('Onset Speed (dist/s)'); xlim([0 11]) %2. onset speed
    subplot(2,3,3); histogram(onsetDist); title('Onset Dist'); xlim([0 2]) %3. onset dist
    subplot(2,3,4); histogram(totalTime); title('Total Time (ms)');xlim([0 1000]) %4. total time
    subplot(2,3,5); histogram(totalSpeed); title('Total Speed (dist/s)'); xlim([0 21]) %5. total speed
    subplot(2,3,6); histogram(totalDist); title('Total Dist'); xlim([0.5 2.5]) %6. total dist
    
    %figure; scatter(onsetTime,onsetDist,10,'filled')
    % Construct matrix for clustering
    attributes = cat(1, onsetTime,onsetSpeed,onsetDist,totalTime,totalSpeed,totalDist)';
    tempAttributes = attributes; tempAttributes(isinf(tempAttributes)) = nan; tempMean = nanmean(tempAttributes,1);
    tempAttributes(badFlag,:) = repmat(tempMean + randi(size(tempMean))*0.001,[sum(badFlag) 1]); 
    infFlag = sum(isinf(tempAttributes),2) >0; tempAttributes(infFlag,:) = repmat(tempMean + randi(size(tempMean))*0.001,[sum(infFlag) 1]); 
    nanFlag = sum(isnan(tempAttributes),2) >0; tempAttributes(nanFlag,:) = repmat(tempMean + randi(size(tempMean))*0.001,[sum(nanFlag) 1]); 
    
    save(['wheelAttributes_' mice{mouse} '.mat'],'tempAttributes');
end 


end



function [clusterProportion] = plotCluster(wheelDownSample_all,RT_all,attribute_tsne,attributes,animalIdx,mice,idx,k)
    figure;
    colors = parula; tempIdx = round(linspace(1, size(colors,1),k+1));
    colors = colors(tempIdx(1:end-1),:);
    plotNTrial = 20;

    % plot the wheel trajectories
    for i = 1:k
        subplot(k,6,(i-1)*6+5)
        tempIdx = find(idx==i); tempIdx = tempIdx(randperm(length(tempIdx))); plotIdx = tempIdx(1:plotNTrial);
        plot(wheelDownSample_all(plotIdx,:)','Color',colors(i,:),'LineWidth',0.2);
        xticks([]);ylim([-1 1]);yticks([-1 0 1]); yticklabels({'R','0','L'});
        tempRT = RT_all(tempIdx);
        title(['Cluster' int2str(i) ' meanRT: ' num2str(nanmean(tempRT),'%.2f') ])
    end
    % plot how learning curve changes with training 
    
    clusterProp = {};  
    for i = 1:length(mice)
        tempIdx = (animalIdx == i); clusterIdx = idx(tempIdx);
        for j = 1:k; clusterProp{i,j} = smoothdata(clusterIdx == j,'movmean',100); end
    end

    clusterProportion = {};
    for i = 1:k
        subplot(k,6,(i-1)*6+6)
        temp = clusterProp(:,i); temp = fn_cell2matFillNan(temp);
        hold on; plot(temp','Color',matlabColors(1,0.3)); plot(nanmean(temp,1),'Color',matlabColors(1))
        xlim([0 2500]); clusterProportion{i} = temp;
    end
    
    [basis, varExp, proj, ~] = fn_pca(attribute_tsne');
    

    subplot(2,6,1); 
    sampleColor = nan(size(attribute_tsne,1),3);
    for i = 1:k; sampleColor(idx==i,:) = repmat(colors(i,:),sum(idx==i),1); end 
    scatter3(proj(1,:),proj(2,:),proj(3,:),4,sampleColor,'filled'); 
    title('pca');hold on;
    
    
    subplot(2,6,7); hold on;
    for i = 1:k; h = cdfplot(RT_all(idx==i));set( h,'Color', colors(i,:)); end
    xlabel('Reaction Time'); ylabel('Cumulative Frequency'); title('Reaction Time'); grid off

    tempPlotList = [2 3 4 8 9 10];
    titleName = {'Onset Time (ms)','Onset Speed (dist/s)','Onset Dist','Total Time (ms)','Total Speed (dist/s)','Total Dist'};
    for j = 1:length(tempPlotList)
        subplot(2,6,tempPlotList(j)); hold on;
        for i = 1:k; h = cdfplot(attributes(idx==i,j));set( h,'Color', colors(i,:)); end
        xlabel(titleName{j}); ylabel('Cumulative Frequency'); title(titleName{j}); grid off;
    end
    
    figure; 
    sampleColor = nan(size(attribute_tsne,1),3);
    for i = 1:k; sampleColor(idx==i,:) = repmat(colors(i,:),sum(idx==i),1); end 
    scatter(proj(1,:),proj(3,:),1,sampleColor,'filled'); 
    title('pca');hold on; xlabel('pc1');  ylabel('pc3')

    figure; 
    subplot(4,1,1);hold on;
    for i = 1:k; h = cdfplot(RT_all(idx==i));set( h,'Color', colors(i,:)); end
    xlabel('Reaction Time (s)'); ylabel('Cumulative Frequency'); grid off
    xlim([0 2.5]); xticks(0:0.5:2.5)

    subplot(4,1,2); hold on;
    for i = 1:k; h = cdfplot(attributes(idx==i,1));set( h,'Color', colors(i,:)); end
    xlabel('Onset Time (ms)'); ylabel('Cumulative Frequency');  grid off;
    xlim([0 2500]); xticks(0:500:2500)

    subplot(4,1,3); hold on;
    for i = 1:k; h = cdfplot(attributes(idx==i,6));set( h,'Color', colors(i,:)); end
    xlim([0 8]); xticks(0:2:8)
    xlabel('Normalized total dist'); ylabel('Cumulative Frequency');  grid off;

    subplot(4,1,4); hold on;
    for i = 1:k; h = cdfplot(attributes(idx==i,6));set( h,'Color', colors(i,:)); end
    xlim([0 8]); xticks(0:2:8)
    xlabel('Normalized total dist'); ylabel('Cumulative Frequency');  grid off;


    


end