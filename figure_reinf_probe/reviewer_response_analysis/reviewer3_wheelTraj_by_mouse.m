%% LOAD DATA
clear; 
mice ={'zz107','zz109','zz111','zz112','zz113','zz115'};
opsParam.biasBlockType = 'AUC';
allMouse = fn_getObjPT_bin40(mice,opsParam); 
mouseMega = wheel2AFCmega(allMouse);
probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
[label,attributesOrig] = fn_loadWheelCluster();
labelCell = {};attributesCell= {};tempTrialSum = [0 cumsum(mouseMega.nTrials)];
for i = 1:mouseMega.nMouse
    labelCell{i} = label(tempTrialSum(i)+1:tempTrialSum(i+1))'; 
    attributesCell{i} = attributesOrig(tempTrialSum(i)+1:tempTrialSum(i+1),:);
    allMouse{i}.behav.label = labelCell{i};
end
mouseMega = wheel2AFCmega(allMouse);

%% Cluster type by animal 
labelCount = nan(4,mouseMega.nMouse);
for i = 1:mouseMega.nMouse
    for j = 1:4
        labelCount(j,i) = sum(labelCell{i}==j);
    end
end

y1 = labelCount./repmat(sum(labelCount,1),size(labelCount,1),1); [~,sortIdx1] = sort(y1(1,:),'descend');

figure; subplot(1,2,1); bar(y1(:,sortIdx1)','stacked'); title('proportion of blocks in first half of trials')

labelCellMat = fn_cell2matFillNan(labelCell);

subplot(1,2,2); hold on;
for i = 1:4
    tempValue = smoothdata(labelCellMat==i,'movmean',100,'includenan');
    tempValue = tempValue ./ repmat(nanmean(smoothdata(~isnan(labelCellMat),'movmean',100),2),[1 size(tempValue,2)]);
    fn_plotMeanErrorbar(1:size(labelCellMat,1), tempValue',matlabColors(i),matlabColors(i),{},{'faceAlpha',0.2});

end
xlim([1 2500])