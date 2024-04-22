clear; 
mice ={'zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);

%% PART 1 --CLUSTERING

%% PART 1.1 -- GET PARAMETERS FOR CLUSTERING

mouse = 6;
[wheelPos, badFlag] = removeBadTrials(allMouse{mouse}.behavVar.wheelPos_aligned);
nTrials = size(wheelPos,1);
behav = allMouse{mouse}.behav;

% normalize the max position of the wheel to 1; downsample every 50ms 
downSampleRate = 50; tempTimeDownSampleIdx = 2001:downSampleRate:4501;
tempMax = wheelPos(:,tempTimeDownSampleIdx(1):tempTimeDownSampleIdx(end)); tempMax = max(tempMax(:));
wheelPos = wheelPos /tempMax;
tempWheelDownSample = wheelPos(:,tempTimeDownSampleIdx);

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
    tempOnsetIdx(i) = find(abs(tempWheelDownSample(i,:))~=0,1); 
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

% save the data of all animals
%tempAttributes = []; animalIdx = []; for i = 1:6; a = load(['wheelAttributes_' mice{i} '.mat']);
%animalIdx( length(animalIdx) + (1:length(a.tempAttributes)) )  = i;
%tempAttributes = cat(1,tempAttributes,a.tempAttributes);
%end
%save('wheelAttributes_allMouse.mat','tempAttributes','animalIdx');
%% PART 1.2 -- DO DIMENSIONALITY REDUCTION FOR VISUALIZATION
%Y = tsne(tempAttributes);
[basis, varExp, proj, covMat] = fn_pca(tempAttributes');
%% PART 1.3 -- DO K-MEANS CLUSTERING
%k = 5; idx = kmeans(attributes,k);
%idx = spectralcluster(attributes,k);
%idx = clusterdata(attributes,'Linkage','centroid','Maxclust',k);
load(['wheel_cluster_' mice{mouse} '.mat'],'attributesScaled','label4')
k=4; idx = label4'+1;
Y = tsne(attributesScaled);
%% PART 1.4 -- VISUALIZE K-MEANS RESULTS
RT = behav.reactionTime;
figure;
colors = parula; tempIdx = round(linspace(1, size(colors,1),k+1));
colors = colors(tempIdx(1:end-1),:);
tempSubplotList = [1 2 4 5 7 8 10 11 13 14];
for i = 1:k
    subplot(5,3,tempSubplotList(i));
    tempIdx = find(idx==i);
    plot(tempWheelDownSample(tempIdx,:)','Color',colors(i,:),'LineWidth',0.2);
    xticks([]);ylim([-1 1]);yticks([-1 0 1]); yticklabels({'R','0','L'});
    tempRT = RT(tempIdx);
    title(['Cluster' int2str(i) ' meanRT: ' num2str(nanmean(tempRT),'%.2f') ])
end

sampleColor = nan(size(attributes,1),3);
for i = 1:k; sampleColor(idx==i,:) = repmat(colors(i,:),sum(idx==i),1); end 
subplot(2,3,3); %scatter(proj(1,:),proj(2,:),30,sampleColor,'filled'); 
scatter(Y(:,1),Y(:,2),30,sampleColor,'filled');
xticks([]); yticks([]); title('t-sne')
% HISTOGRAM OF K-MEANS CLUSTERS
subplot(2,3,6); hold on;
for i = 1:k; h = cdfplot(RT(idx==i));set( h,'Color', colors(i,:)); end
xlabel('Reaction Time'); ylabel('Cumulative Frequency'); title('Reaction Time')

% CLUSTER ACROSS LEARNING
figure;
for i = 1:10
   subplot(5,2,i);
   plot(smoothdata(idx==i,'movmean',300));
   tempLRT = nanmean(behav.reactionTime(idx==i & behav.bias>0.2)); tempRRT = nanmean(behav.reactionTime(idx==i & behav.bias<-0.2));
   tempLperc = sum(idx==i & behav.bias>0.2) / sum(idx==i); tempRperc = sum(idx==i & behav.bias<-0.2) / sum(idx==i);
   tempAcc = nanmean(behav.responseType(idx==i)==1);
   title(['RT ' num2str(tempLRT,'%.2f') ' ' num2str(tempRRT,'%.2f') ...
       ' perc ' num2str(tempLperc,'%.2f') ' ' num2str(tempRperc,'%.2f') ' acc ' num2str(tempAcc,'%.2f') ])
end



%% PART 2 -- REINF VS. PROBE COMPARISON
reinfTrajL = tempWheelDownSample((~behav.goodProbe) & behav.action==1,:);
reinfTrajR = tempWheelDownSample((~behav.goodProbe) & behav.action==2,:);

figure;subplot(2,1,1);  hold on;
plot(0:50,tempWheelDownSample(behav.goodProbe & behav.action==1,:),'Color',matlabColors(1,0.4));
plot(0:50,nanmean(tempWheelDownSample(behav.goodProbe & behav.action==1,:),1),'Color',matlabColors(1),'Linewidth',3);
plot(0:50,nanmean(reinfTrajL,1),'Color',[0.6 0.6 0.6],'Linewidth',3); xlim([0 50]);
title(['mouse ' mice{mouse}])
subplot(2,1,2);hold on; 
plot(0:50,tempWheelDownSample(behav.goodProbe & behav.action==2,:),'Color',matlabColors(1,0.4));
plot(0:50,nanmean(tempWheelDownSample(behav.goodProbe & behav.action==2,:),1),'Color',matlabColors(1),'Linewidth',3);
plot(0:50,nanmean(reinfTrajR,1),'Color',[0.6 0.6 0.6],'Linewidth',3); xlim([0 50]);

% look at probe vs. reinf only at the points with significant difference
probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
outCell = mouseMega.objFun('binProbeByTrialFromLearningOnset',{[reinfThre probeThre],probeTrialBin});
 
reinfAlignPoint = cell2mat(outCell{3});
reinfSelStart = max(reinfAlignPoint(mouse)-600,1);
reinfSelEnd = reinfAlignPoint(mouse)+600;

tempBehav = behav(reinfSelStart:reinfSelEnd,:);
tempWheelDownSampleSel = tempWheelDownSample(reinfSelStart:reinfSelEnd,:);
figure;subplot(2,1,1);  hold on;
plot(0:50,tempWheelDownSampleSel(tempBehav.goodProbe & tempBehav.action==1,:),'Color',matlabColors(1,0.4));
plot(0:50,nanmean(tempWheelDownSampleSel(tempBehav.goodProbe & tempBehav.action==1,:),1),'Color',matlabColors(1),'Linewidth',3);
plot(0:50,nanmean(reinfTrajL,1),'Color',[0.6 0.6 0.6],'Linewidth',3); xlim([0 50]);
title(['mouse ' mice{mouse}])
subplot(2,1,2);hold on; 
plot(0:50,tempWheelDownSampleSel(tempBehav.goodProbe & tempBehav.action==2,:),'Color',matlabColors(1,0.4));
plot(0:50,nanmean(tempWheelDownSampleSel(tempBehav.goodProbe & tempBehav.action==2,:),1),'Color',matlabColors(1),'Linewidth',3);
plot(0:50,nanmean(reinfTrajR,1),'Color',[0.6 0.6 0.6],'Linewidth',3); xlim([0 50]);

%% FUNCTIONS
function [tempWheel, badFlag] = removeBadTrials(tempWheel)
    totalTrial = size(tempWheel,1);
    badFlag = tempWheel(:,end)==0;
    tempWheel(badFlag,:) = nan;
    disp(['Bad Trials = ' int2str(sum(badFlag)) ' out of ' int2str(totalTrial)])
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
