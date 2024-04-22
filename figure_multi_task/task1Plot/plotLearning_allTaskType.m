%%
clear; 
mice ={ {'zz054','zz062','zz063','zz066','zz067','zz068','zz069'},...
         {'zz107','zz109','zz111','zz112','zz113','zz115'},...
        {'zz101','zz102','zz105'},...
        {'zz097','zz098'},...
        {'zz121','zz122','zz123'}...
        }; 
nGroup = length(mice);
taskName = {'puretone','puretone new','FM Dir Dur pip','PT pip','FM Dir pip'};
mouseMega = {}; allMouse = {};
for i=1:nGroup
    allMouse{i} = fn_getObjTask1(mice{i});    
end
%%
for i = 1:nGroup
    for j = 1:length(mice{i})
        allMouse{i}{j}.ops.learningCurveBin  = 1200;
        allMouse{i}{j} = allMouse{i}{j}.getAcc();
    end
    mouseMega{i} = wheel2AFCmega(allMouse{i});
end
%% 
% GET PEAK IN THE LEARNING CURVES
peakIdx = {}; peakPerf = {};
figure;
for i = 1:nGroup
    reinf = mouseMega{i}.loadReinf();
    for j = 1: size(reinf.acc,2)
        temp = reinf.acc(:,j); lastidx = find(~isnan(reinf.acc(:,j)),1,'last'); 
        temp(lastidx+1:lastidx+50) = 0.65;
        [pks,p] = islocalmax (temp); 
        subplot(nGroup,size(reinf.acc,2),(i-1)*size(reinf.acc,2)+j)
        flag = p>0.02;
        plot(reinf.acc(:,j)); hold on; 
        scatter((find(flag)), reinf.acc(flag,j),20,'filled');
        
        if sum(flag)==1; peakIdx{i,j} = find(flag); peakPerf{i,j} = reinf.acc(flag,j);
        elseif sum(flag)>1 % if there are more than 1 peak 
            % take the first index that is less than 0.1 accuracy away from
            % the maximum accuracy peak
            tempPerf = reinf.acc(flag,j); tempIdx = find(flag); tempEntry = find(tempPerf > max(tempPerf)*0.9,1,'first');
            peakIdx{i,j} = tempIdx(tempEntry); peakPerf{i,j} = reinf.acc(peakIdx{i,j},j);
        else
            peakIdx{i,j} = nan; peakPerf{i,j} = nan;
        end
        scatter(peakIdx{i,j}, peakPerf{i,j},20,'filled');
        
    end 
end
% GET DROP IN THE LEARNING CURVES
dropIdx = {};dropPerf = {};
figure;
for i = 1:nGroup
    reinf = mouseMega{i}.loadReinf();
    for j = 1: size(reinf.acc,2)
        temp = reinf.acc(:,j); lastidx = find(~isnan(reinf.acc(:,j)),1,'last'); 
        temp(lastidx+1:lastidx+50) = 0.75;
        [pks,p] = islocalmin (temp); 
        subplot(nGroup,size(reinf.acc,2),(i-1)*size(reinf.acc,2)+j)
        flag = p>0.02;
        plot(reinf.acc(:,j)); hold on; 
        scatter((find(flag)), reinf.acc(flag,j),20,'filled');
        
        if sum(flag)==1; tempIdx = find(flag); 
            if tempIdx > peakIdx{i,j} && peakPerf{i,j} - reinf.acc(flag,j) > 0.15
                dropIdx{i,j} = tempIdx; dropPerf{i,j} = reinf.acc(flag,j); 
            else
                dropIdx{i,j} = nan; dropPerf{i,j} = nan;
            end
        elseif sum(flag)>1 % if there are more than 1 peak 
            % take the first index that is less than 0.1 accuracy away from
            % the maximum accuracy peak
            tempIdx = find(flag); tempEntry = find(tempIdx > peakIdx{i,j});
            tempAccDrop = repmat(peakPerf{i,j},[length(tempEntry) 1]) - reinf.acc(tempIdx(tempEntry),j);
            tempSelect = find(tempAccDrop>0.15,1,'first'); tempIdx = tempIdx(tempEntry(tempSelect));
            dropIdx{i,j} = tempIdx; dropPerf{i,j} = reinf.acc(dropIdx{i,j},j);
        else       
            dropIdx{i,j} = nan; dropPerf{i,j} = nan;
        end
        if isempty(dropIdx{i,j});dropIdx{i,j} = nan; dropPerf{i,j} = nan; end
        if ~isempty(dropIdx{i,j}); scatter(dropIdx{i,j}, dropPerf{i,j},20,'filled'); end 
    end 
end

%% Look at the difference of learning curve
flatIdx = {}; decayIdx = {}; flatPerf = {}; decayPerf= {};
figure;
for i = 1:nGroup
    reinf = mouseMega{i}.loadReinf();
    for j = 1: size(reinf.acc,2)
        subplot(nGroup,size(reinf.acc,2),(i-1)*size(reinf.acc,2)+j)
        temp = reinf.acc(:,j); temp = smoothdata(temp,'movmean',400,'includenan'); tempDiff = diff(temp); tempDiff = smoothdata(tempDiff,'movmean',400,'includenan');
        flatIdx{i,j} =  find(tempDiff < 0.5e-4 & reinf.acc(1:end-1,j)>0.75); flatIdx{i,j} =  flatIdx{i,j}(find(flatIdx{i,j} < peakIdx{i,j},1,'first')); 
        
        decayIdx{i,j} =  find(tempDiff < -1e-4); decayIdx{i,j} =  decayIdx{i,j}(find(decayIdx{i,j} > peakIdx{i,j},1,'first')); 
        
        if isempty(flatIdx{i,j}); flatIdx{i,j} = peakIdx{i,j}; end 
        plot(tempDiff); hold on;
        tempFirstIdx = find(~isnan(tempDiff),1,'first'); tempLastIdx = find(~isnan(tempDiff),1,'last');
        plot([tempFirstIdx tempLastIdx],[0 0],'Color',[0.8 0.8 0.8]); xlim([tempFirstIdx tempLastIdx]);
        flatPerf{i,j} = reinf.acc(flatIdx{i,j},j); decayPerf{i,j} = reinf.acc(decayIdx{i,j},j);
        if isempty(decayIdx{i,j}); decayIdx{i,j} = nan; end 
        if isempty(decayPerf{i,j}); decayPerf{i,j} = nan; end 
    end
end


%%
figure; 
for i = 1:nGroup
    reinf = mouseMega{i}.loadReinf();
    for j = 1: size(reinf.acc,2)
        subplot(nGroup,size(reinf.acc,2),(i-1)*size(reinf.acc,2)+j)
        plot(reinf.acc(:,j)); hold on;
        
        scatter(flatIdx{i,j}, flatPerf{i,j},20,matlabColors(2),'filled'); 
        scatter(decayIdx{i,j}, decayPerf{i,j},20,matlabColors(2),'filled'); 
        %scatter(peakIdx{i,j}, peakPerf{i,j},20,matlabColors(4),'filled'); 
        scatter(dropIdx{i,j}, dropPerf{i,j},20,matlabColors(5),'filled'); 
    end
end


%% plot the peak trial index
figure; subplot(1,2,1); plotBarTasks(peakIdx,taskName);xtickangle(45); title('Trials to Peak')
subplot(1,2,2); plotBarTasks(peakPerf,taskName);xtickangle(45); title('Accuracy at Peak')

figure; subplot(2,3,1); plotBarTasks(flatIdx,taskName); xtickangle(45); title('Trials to Plateau')
subplot(2,3,2); plotBarTasks(flatPerf,taskName);xtickangle(45); title('Accuracy at Plateau')
subplot(2,3,3); scatter(cell2mat(flatIdx(:)),cell2mat(flatPerf(:)),'filled'); title('Trials vs. performance correlation')

subplot(2,3,4);plotBarTasksDiff(decayIdx,flatIdx,taskName);xtickangle(45); title('Trials at Plateau')
subplot(2,3,5);plotBarTasksDiff(peakPerf,dropPerf,taskName);xtickangle(45); title('Performance drop after plateau')
subplot(2,3,6); scatter(cell2mat(decayIdx(:))-cell2mat(flatIdx(:)),cell2mat(peakPerf(:))-cell2mat(dropPerf(:)),'filled'); 
title('Trials at plateau vs. performance drop')
%% Use an arbitary threshold for 
criteriaOn = {}; criteriaOff = {}; criteriaLen = {}; accuracyThre = 0.75;
for i = 1:nGroup
    reinf = mouseMega{i}.loadReinf();
    for j = 1: size(reinf.acc,2)
        temp = reinf.acc(:,j); 
        [startIdx,endIdx,len] = fn_findBlock(temp > accuracyThre,1,50);
        if ~isempty(startIdx); criteriaOn{i,j} = startIdx(1); criteriaOff{i,j} = endIdx(1); criteriaLen{i,j} = len(1);
        else;  criteriaOn{i,j} = nan; criteriaOff{i,j} = nan; criteriaLen{i,j} = nan; end
        if find(~isnan(temp),1,'Last') == criteriaOff{i,j}; criteriaOff{i,j} = nan; end 
    end
end

figure; subplot(2,1,1); plotBarTasks(criteriaOn,taskName); xtickangle(45);title(['Trials to criteria=' num2str(accuracyThre,'%.2f')])
subplot(2,1,2); plotBarTasksDiff(criteriaOff,criteriaOn,taskName); xtickangle(45);title(['Trials above criteria at plateau=' num2str(accuracyThre,'%.2f')])
%% plot the learning curve of all animals

for i = 1:nGroup
    for j = 1:length(mice{i})
        allMouse{i}{j}.ops.learningCurveBin  = 100;
        allMouse{i}{j} = allMouse{i}{j}.getAcc();
    end
    mouseMega{i} = wheel2AFCmega(allMouse{i});
end

figure; trialLim = [1 6000]; tempAcc = nan(trialLim(end),nGroup);
for i = 1:nGroup
    subplot(nGroup+1,1,i); hold on 
    reinf = mouseMega{i}.loadReinf();
    plot(reinf.acc,'Color',[0.8 0.8 0.8]);
    plot(nanmean(reinf.acc,2),'Color',matlabColors(i));
    xlim(trialLim); ylim([0.4 1]); xticks(0:1000:6000)
    if size(reinf.acc,1) < trialLim(end)
        tempAcc(1:size(reinf.acc,1),i) = nanmean(reinf.acc,2);
    else
        tempAcc(:,i) = nanmean(reinf.acc(trialLim(1):trialLim(2),:),2);
    end
end
subplot(nGroup+1,1,nGroup+1);
plot(tempAcc); xlim(trialLim);  xticks(0:1000:6000)

%% Plot functions
function tempMean = plotBarTasks(mat,taskName)
 hold on; tempMean = [];
for i = 1:size(mat,1); tempMean(i) = nanmean([mat{i,:}]); end
bar(tempMean,'EdgeColor',matlabColors(2),'FaceColor',[1 1 1],'LineWidth',2)
for i = 1:size(mat,1);  scatter(i*ones(size([mat{i,:}])), [mat{i,:}],20,[0.4 0.4 0.4],'filled'); end
xticks(1:5); xticklabels(taskName);
end

function tempMean = plotBarTasksDiff(mat1,mat2,taskName)
 hold on; tempMean = [];
for i = 1:size(mat1,1); tempMean(i) = nanmean([mat1{i,:}] - [mat2{i,:}]); end
bar(tempMean,'EdgeColor',matlabColors(2),'FaceColor',[1 1 1],'LineWidth',2)
for i = 1:size(mat1,1);  scatter(i*ones(size([mat1{i,:}])), [mat1{i,:}] - [mat2{i,:}],20,[0.4 0.4 0.4],'filled'); end
xticks(1:5); xticklabels(taskName);
end
