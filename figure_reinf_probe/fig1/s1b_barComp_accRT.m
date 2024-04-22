%% LOAD ANIMALS
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);
reinf = mouseMega.loadReinf;

nMouse = size(reinf.acc,2);
plotLim = 3000;

%% PLOT 1 -- Learning Curve
wheelLim = [35,35,35,35,35,35];
fn_figureSmartDim('hSize',0.4,'widthHeightRatio',0.6); 
tempWheel = getWheel('choiceMoveCorr',mouseMega);
tempWheel = fn_cell2matFillNan(tempWheel);
tempWheel = tempWheel ./ repmat(wheelLim,[size(tempWheel,1) 1]);
tempWheel = smoothdata(tempWheel,'movmean',100,'includenan');
subplot(2,1,1); plotIndividualReinf(tempWheel(1:plotLim,:));
xlim([1 plotLim]); ylabel('correct wheel move'); ylim([-1 4])

subplot(2,1,2);hold on;
plotIndividualReinf(reinf.actionRate(1:plotLim,:));
xlim([1 plotLim]); ylabel('actionRate'); ylim([0 1.0])

%% PLOT 2 -- Naive vs. Expert Bar Comparison

tempBin = 300; start = 1;
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

% get mean acc and RT
tempBin = 300; start = 1; preLearningTrial = start:start+tempBin-1;
actionPre = cellfun(@(x)(mean(x.behav.actionRate(preLearningTrial))),...
    mouseMega.mouseCell);
actionPost = cellfun(@(x,y)(mean(x.behav.actionRate(y))),...
    mouseMega.mouseCell,postLearningTrial);

wheelPre = nan(1,size(tempWheel,2)); wheelPost = nan(1,size(tempWheel,2));
for i = 1:size(tempWheel,2)
    wheelPre(i) = nanmean(tempWheel(preLearningTrial,i));
    wheelPost(i) = nanmean(tempWheel(postLearningTrial{i+7},i));
end
figure;
subplot(2,1,1);
[p] = fn_plotComparison({wheelPre,wheelPost},'paired',true,'dotType','side','compType','errorbarWithDot','scatterArgIn',...
    {5,'MarkerEdgeColor','none','MarkerFaceColor','none'},'errorbarArgIn', {'Color',[0.2 0.2 0.2],'LineWidth',1.5,'LineStyle','none'}); 
xlim([0.6 2.4]); ylim([-10 100]); title(['correct wheel move, p=' num2str(p)])

subplot(2,1,2);
[p] = fn_plotComparison({actionPre,actionPost},'paired',true,'dotType','side','compType','errorbarWithDot','scatterArgIn',...
    {5,'MarkerEdgeColor','none','MarkerFaceColor','none'},'errorbarArgIn', {'Color',[0.2 0.2 0.2],'LineWidth',1.5,'LineStyle','none'}); 
xlim([0.6 2.4]); ylim([0 1.0]); yticks(0:0.5:1.0); title(['actionRate, p=' num2str(p)])

%% all functions
function plotIndividualReinf(mat)
    % PLOT REINF
    hold on;
    plot(mat,'Color', fn_wheelColorsPT('Reinf',0.3), 'LineWidth',1);
    plot(nanmean(mat,2),'Color', fn_wheelColorsPT('Reinf'), 'LineWidth',3);
end


function tempWheel = getWheel(attrName,mouseMega)       
    tempWheel = {};
    for i = 1:mouseMega.nMouse
        if ismember(attrName, mouseMega.mouseCell{i}.behav.Properties.VariableNames)
            temp = mouseMega.mouseCell{i}.behav.(attrName);
            tempWheel = [tempWheel; {temp}];
        end 
    end 
end 