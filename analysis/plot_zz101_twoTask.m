%% GET OBJECT FOR PT pip COHORTS
clear; 
mouse ='zz101';
task = 'Two_task';

rootPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\';
loadPath = [rootPath filesep 'trialData\'];
load([loadPath mouse '.mat']);

mouseObj = wheel2AFC(trialData,'opsPath',['C:\Users\zzhu34\Documents\tempdata\octoData\opsFiles\' task ],'mouse',mouse );

%mouseObj = mouseObj.removeMiss(); 

mouseObj = mouseObj.getAcc();


%% GET DATA 

task_names = {'FM_Dir_Dur_pip','PT_pip','Two_task'};
taskIdx= {}; taskDay = {};
for i = 1:length(task_names)
    tempIdx = strcmp(mouseObj.behav.trainingType, task_names{i});
    taskIdx{i} = tempIdx;
    taskDay{i} = unique(mouseObj.behav.day(tempIdx));
end

% EXPERT LEVEL ACCURACY FOR THE TWO TASKS 
trialBin = 80;
taskAcc = {}; taskBias = {}; taskRT = {};
for i = 1:2 % loop through both tasks
    % take the accuracy of the last two days
    dayIdx = mouseObj.behav.day == taskDay{i}(end);%| mouseObj.behav.day == taskDay{i}(end-1);
    tempIdx = taskIdx{i} & dayIdx;
    tempStim = mouseObj.behav.stimulus(tempIdx);
    tempResp = mouseObj.behav.responseType(tempIdx);
    
    [bias,acc_L,acc_R] = fn_getBias(tempStim,tempResp,trialBin);

    if i == 1
        acc_L = acc_L(150:end); acc_R = acc_R(150:end);
        bias = bias(150:end);
    end
    taskAcc{i,1} = nanmean((acc_L+acc_R) /2); taskBias{i,1} = nanmean(abs(bias));
    taskAcc{i,2} = nanstd((acc_L+acc_R)/2); % std
    taskBias{i,2} = nanstd(abs(bias)); % std
end

% GET INTERLEAVED TRAINING DATA

[behav, task1Idx,task2Idx,dayIdx] = selectInterLeaveDay(mouseObj,taskDay{3});

taskAcc1 = smoothdata(behav.responseType(task1Idx)==1,'movmean',trialBin)./...
    smoothdata(behav.responseType(task1Idx)~=0,'movmean',trialBin); 
taskAcc2 = smoothdata(behav.responseType(task2Idx)==1,'movmean',trialBin)./...
    smoothdata(behav.responseType(task1Idx)~=0,'movmean',trialBin); 


% GET
[ctxtFlag, ctxtIdx,blockBin,blockCtxt,blockDay]=getContextIdx(dayIdx);
taskAcc1_block = []; taskAcc2_block = [];taskBias1_block = [];taskBias2_block = [];
taskRT1_block = []; taskRT2_block = [];
for i = 1:length(blockBin)-1
    trialIdx = (blockBin(i)+1):blockBin(i+1);
    switch blockCtxt(i)
        case 1
            [bias,acc] = fn_getBiasMean(behav.stimulus(trialIdx),behav.responseType(trialIdx));
            taskAcc1_block(i) = acc; taskAcc2_block(i) = nan; 
            taskBias1_block(i) = abs(bias);  taskBias2_block(i) = nan; 
            taskRT1_block(i) = mean(behav.reactionTime(trialIdx)); taskRT2_block(i) = nan; 
        case 2
            [bias,acc] = fn_getBiasMean(behav.stimulus(trialIdx)-2,behav.responseType(trialIdx));
            taskAcc1_block(i) = nan; taskAcc2_block(i) = acc; 
            taskBias1_block(i) = nan;  taskBias2_block(i) = abs(bias); 
            taskRT1_block(i) = nan; taskRT2_block(i) = mean(behav.reactionTime(trialIdx)); 
        case 3
            tempStim = behav.stimulus(trialIdx); tempResp = behav.responseType(trialIdx);
            tempRT = behav.reactionTime(trialIdx);
            task1Flag = tempStim==1 | tempStim==2; task2Flag = tempStim==3 | tempStim==4;

            [bias,acc] = fn_getBiasMean(tempStim(task1Flag),tempResp(task1Flag));
            taskAcc1_block(i) = acc;  taskBias1_block(i) = abs(bias);  
            taskRT1_block(i) = mean(tempRT(task1Flag));

            [bias,acc] = fn_getBiasMean(tempStim(task2Flag)-2,tempResp(task2Flag));
            taskAcc2_block(i) = acc; taskBias2_block(i) = abs(bias);
            taskRT2_block(i) = mean(tempRT(task2Flag));
    end
end
taskRT1_block(134) = nan; taskRT2_block(134) = nan;
taskRT1_block(73:86) = nan; taskRT2_block(73:86) = nan;


%% PLOT 1.1 -- PLOT PERFORMANCE OF TWO TASKS WITH EXPERT LEVEL PERFORMANCE
figure; hold on; 
xlimit = max([length(taskAcc1), length(taskAcc2)]);
ylimit = [floor(min(taskAcc1)*10)/10 1];
plot(taskAcc1,'LineWidth',2);  plot(taskAcc2,'LineWidth',2); 

tempX = [1 xlimit xlimit 1];
fn_plotFillErrorbar([1 xlimit],[taskAcc{1,1} taskAcc{1,1}],[taskAcc{1,2} taskAcc{1,2}],matlabColors(1,0.8),...
    'faceAlpha',0.2,'LineStyle','none');
fn_plotFillErrorbar([1 xlimit],[taskAcc{2,1} taskAcc{2,1}],[taskAcc{2,2} taskAcc{2,2}],matlabColors(2,0.8),...
    'faceAlpha',0.2,'LineStyle','none');

for i = 1:length(dayIdx)-1
    plot([round(dayIdx(i)/2) round(dayIdx(i)/2)],ylimit,'Color',[0.8 0.8 0.8]);
end

ylim(ylimit);xlim([1 xlimit])

%% PLOT 1.2 -- PLOT PERFORMANCE OF TWO TASKS WITH COLORED BLOCKS (NOT ACCURATE)

figure; hold on; 
xlimit = max([length(taskAcc1), length(taskAcc2)]);
ylimit = [floor(min(taskAcc1)*10)/10 1];
plot(taskAcc1,'LineWidth',2);  plot(taskAcc2,'LineWidth',2); 

for i = 1:length(ctxtIdx)% loop through context
    for j = 1:size(ctxtIdx{i},1)
        fn_plotFillErrorbar(round([ctxtIdx{i}(j,1) ctxtIdx{i}(j,2)]/2),[mean(ylimit) mean(ylimit)],...
            [(ylimit(2)-ylimit(1))/2 (ylimit(2)-ylimit(1))/2],matlabColors(2+i,0.8),...
            'faceAlpha',0.2,'LineStyle','none');       
    end     
end
ylim(ylimit);xlim([1 xlimit])


%% PLOT 1.3 -- PLOT PERFORMANCE OF FIRST FEW DAYS FO THE TRANSFER
%selectDay = [taskDay{1}(end-1:end);taskDay{3}(1:end)];
selectDay = [taskDay{3}(1:13)];
[behav, task1Idx,task2Idx,dayTransIdx] = selectInterLeaveDay(mouseObj,selectDay);
taskAcc1 = {}; taskAcc2 = {}; taskBias1 = {}; taskBias2 = {}; taskRT1 = {}; taskRT2 = {};
for i = 1:length(selectDay)

    tempIdx  = behav.day==selectDay(i);
    taskAcc1{i} = smoothdata(behav.responseType(tempIdx & task1Idx)==1,'movmean',trialBin)./...
        smoothdata(behav.responseType(tempIdx & task1Idx)~=0,'movmean',trialBin); 
    taskAcc2{i} = smoothdata(behav.responseType(tempIdx & task2Idx)==1,'movmean',trialBin)./...
        smoothdata(behav.responseType(tempIdx & task2Idx)~=0,'movmean',trialBin); 
    
    taskBias1{i} = fn_getBias( behav.stimulus(tempIdx & task1Idx), behav.responseType(tempIdx & task1Idx), trialBin); 
    taskBias2{i} = fn_getBias( behav.stimulus(tempIdx & task2Idx)-2, behav.responseType(tempIdx & task2Idx), trialBin); 

    taskRT1{i} = smoothdata(behav.reactionTime(tempIdx & task1Idx),'movmean',trialBin);
    taskRT2{i} = smoothdata(behav.reactionTime(tempIdx & task2Idx),'movmean',trialBin); 
end
%taskAcc2{1} = nan(size(taskAcc1{1})); taskAcc2{2} = nan(size(taskAcc1{2}));
%taskBias2{1} = nan(size(taskBias1{1})); taskBias2{2} = nan(size(taskBias1{2}));
%taskRT2{1} = nan(size(taskRT1{1})); taskRT2{2} = nan(size(taskRT1{2}));

% plot accuracy
figure; hold on; ylimit = [0.3 1]; tempCount = 0;  
for i = 1:length(selectDay)
    plot((tempCount+1):(tempCount+length(taskAcc1{i})), taskAcc1{i}, 'Color',matlabColors(1),'LineWidth',2); 
    plot((tempCount+1):(tempCount+length(taskAcc1{i})),taskAcc2{i}, 'Color',matlabColors(2),'LineWidth',2);
    tempCount = tempCount+length(taskAcc1{i});
    plot([tempCount tempCount],ylimit,'Color',[0.8 0.8 0.8]);
end
plot([1 tempCount],[0.5 0.5],'Color',[0.8 0.8 0.8]);
fn_plotFillErrorbar([1 tempCount],[taskAcc{1,1} taskAcc{1,1}],[taskAcc{1,2} taskAcc{1,2}],matlabColors(1,0.8),...
    'faceAlpha',0.2,'LineStyle','none');
fn_plotFillErrorbar([1 tempCount],[taskAcc{2,1} taskAcc{2,1}],[taskAcc{2,2} taskAcc{2,2}],matlabColors(2,0.8),...
    'faceAlpha',0.2,'LineStyle','none');
ylim(ylimit);xlim([1 tempCount])

% plot RT
figure; hold on; ylimit = [0 2.5]; tempCount = 0;  
for i = 1:length(selectDay)
    plot((tempCount+1):(tempCount+length(taskRT1{i})), taskRT1{i}, 'Color',matlabColors(1),'LineWidth',2); 
    plot((tempCount+1):(tempCount+length(taskRT1{i})),taskRT2{i}, 'Color',matlabColors(2),'LineWidth',2);
    tempCount = tempCount+length(taskRT1{i});
    plot([tempCount tempCount],ylimit,'Color',[0.8 0.8 0.8]);
end
plot([1 tempCount],[0.5 0.5],'Color',[0.8 0.8 0.8]);
ylim(ylimit);xlim([1 tempCount])

% plot bias
figure; hold on; ylimit = [-1 1]; tempCount = 0;  
for i = 1:length(selectDay)
    plot((tempCount+1):(tempCount+length(taskBias1{i})), taskBias1{i}, 'Color',matlabColors(1),'LineWidth',2); 
    plot((tempCount+1):(tempCount+length(taskBias2{i})),taskBias2{i}, 'Color',matlabColors(2),'LineWidth',2);
    tempCount = tempCount+length(taskBias1{i});
    plot([tempCount tempCount],ylimit,'Color',[0.8 0.8 0.8]);
end
plot([1 tempCount],[0 0],'Color',[0.8 0.8 0.8]);
fn_plotFillErrorbar([1 tempCount],[taskBias{1,1} taskBias{1,1}],[taskBias{1,2} taskBias{1,2}],matlabColors(1,0.8),...
    'faceAlpha',0.2,'LineStyle','none');
fn_plotFillErrorbar([1 tempCount],[taskBias{2,1} taskBias{2,1}],[taskBias{2,2} taskBias{2,2}],matlabColors(2,0.8),...
    'faceAlpha',0.2,'LineStyle','none');
ylim(ylimit);xlim([1 tempCount])


% plot bias
taskBias1 = fn_getBias( behav.stimulus(task1Idx), behav.responseType(task1Idx), trialBin); 
taskBias2 = fn_getBias( behav.stimulus(task2Idx)-2, behav.responseType(task2Idx), trialBin); 
%taskBias2 = cat(1,nan(600,1),taskBias2);

figure; hold on; ylimit = [-1 1]; 
%tempIdx = diff([0 dayTransIdx]); tempIdx(3:end) = tempIdx(3:end)/2; tempIdx = cumsum(tempIdx);
tempIdx = diff([0 dayTransIdx]); tempIdx = tempIdx/2; tempIdx = cumsum(tempIdx);
for i = 1:length(tempIdx)
    %plot((tempCount+1):(tempCount+length(taskBias1{i})), taskBias1{i}, 'Color',matlabColors(1),'LineWidth',2); 
    %plot((tempCount+1):(tempCount+length(taskBias2{i})),taskBias2{i}, 'Color',matlabColors(2),'LineWidth',2);
    plot([tempIdx(i) tempIdx(i)],ylimit,'Color',[0.8 0.8 0.8]);
end
plot(taskBias1, 'Color',matlabColors(1),'LineWidth',2);plot(taskBias2, 'Color',matlabColors(2),'LineWidth',2);
plot([1 tempIdx(end)],[0 0],'Color',[0.8 0.8 0.8]);
fn_plotFillErrorbar([1 tempCount],[taskBias{1,1} taskBias{1,1}],[taskBias{1,2} taskBias{1,2}],matlabColors(1,0.8),...
    'faceAlpha',0.2,'LineStyle','none');
fn_plotFillErrorbar([1 tempIdx(end)],[taskBias{2,1} taskBias{2,1}],[taskBias{2,2} taskBias{2,2}],matlabColors(2,0.8),...
    'faceAlpha',0.2,'LineStyle','none');
ylim(ylimit); xlim([1 length(taskBias1)])



%% PLOT 2.1 -- PLOT THE REACTION TIME 
xdata = 1:length(blockCtxt);
figure; hold on;
tempIsnan = isnan(taskAcc1_block); plot(xdata(~tempIsnan), taskAcc1_block(~tempIsnan))
tempIsnan = isnan(taskAcc2_block); plot(xdata(~tempIsnan), taskAcc2_block(~tempIsnan))


figure; hold on;
tempIsnan = isnan(taskRT1_block); plot(xdata(~tempIsnan), taskRT1_block(~tempIsnan))
tempIsnan = isnan(taskRT2_block); plot(xdata(~tempIsnan), taskRT2_block(~tempIsnan))


figure; subplot(1,2,1); hold on;
plot(ones(1,sum(blockCtxt~=3)), taskRT1_block(blockCtxt~=3),'o')
plot(ones(1,sum(blockCtxt==3))*2, taskRT1_block(blockCtxt==3),'o')

subplot(1,2,2); hold on
plot(ones(1,sum(blockCtxt~=3)), taskRT2_block(blockCtxt~=3),'o')
plot(ones(1,sum(blockCtxt==3))*2, taskRT2_block(blockCtxt==3),'o')

% select out the first and last block of each day
startDayFlag = logical([1 diff(blockDay)]); endDayFlag = logical([diff(blockDay) 1]);
taskRT1_block (startDayFlag | endDayFlag) = nan;
taskRT2_block (startDayFlag | endDayFlag) = nan;

figure; subplot(1,2,1); hold on;
plot(ones(1,sum(blockCtxt~=3)), taskRT1_block(blockCtxt~=3),'o')
plot(ones(1,sum(blockCtxt==3))*2, taskRT1_block(blockCtxt==3),'o')

subplot(1,2,2); hold on
plot(ones(1,sum(blockCtxt~=3)), taskRT2_block(blockCtxt~=3),'o')
plot(ones(1,sum(blockCtxt==3))*2, taskRT2_block(blockCtxt==3),'o')



taskAcc1_block (startDayFlag | endDayFlag) = nan;
taskAcc2_block (startDayFlag | endDayFlag) = nan;
%% combine blocks within a day

selectDay = [taskDay{3}(1:13)];
[behav, task1Idx,task2Idx,dayIdx] = selectInterLeaveDay(mouseObj,selectDay);

tempDayIdx = [0 dayIdx]; rmTrials = 80;

for i = 1:length(tempDayIdx)-1
    tempIdx = (tempDayIdx(i)+1+rmTrials):(tempDayIdx(i+1)-rmTrials);
    tempRT = behav.reactionTime(tempIdx); 
    tempStim = behav.stimulus(tempIdx)'; 

    block1 = ctxtFlag(tempIdx) == 1; taskRT1_day(i,1) = mean(tempRT(block1));
    block2 = ctxtFlag(tempIdx) == 2; taskRT2_day(i,1) = mean(tempRT(block2));
    block3 = ctxtFlag(tempIdx) == 3; 
    taskRT1_day(i,2) = mean(tempRT(block3 & (tempStim==1 | tempStim==2)));
    taskRT2_day(i,2) = mean(tempRT(block3 & (tempStim==3 | tempStim==4)));

end

figure; subplot(1,2,1); plot(taskRT1_day)
subplot(1,2,2); plot(taskRT2_day)


figure; subplot(1,2,1); plot(taskRT1_day(:,1) -taskRT1_day(:,2))
subplot(1,2,2); plot(taskRT2_day(:,1) - taskRT2_day(:,2))

%% combine acc within a day

%selectDay = [taskDay{3}(17:18);taskDay{3}(20:22)];
selectDay = [taskDay{3}(1:end)];
[behav, task1Idx,task2Idx,dayIdx] = selectInterLeaveDay(mouseObj,selectDay);

tempDayIdx = [0 dayIdx]; rmTrials = 0;

taskAcc1_day = []; taskAcc2_day = [];

for i = 1:length(tempDayIdx)-1
    tempIdx = (tempDayIdx(i)+1+rmTrials):(tempDayIdx(i+1)-rmTrials);
    tempResp = behav.responseType(tempIdx); 
    tempStim = behav.stimulus(tempIdx)'; 

    block1 = ctxtFlag(tempIdx) == 1; taskAcc1_day(i,1) = mean(tempResp(block1)==1);
    block2 = ctxtFlag(tempIdx) == 2; taskAcc2_day(i,1) = mean(tempResp(block2)==1);
    block3 = ctxtFlag(tempIdx) == 3; 
    taskAcc1_day(i,2) = mean(tempResp(block3 & (tempStim==1 | tempStim==2))==1);
    taskAcc2_day(i,2) = mean(tempResp(block3 & (tempStim==3 | tempStim==4))==1);

end

figure; subplot(1,3,1); hold on; temp = taskAcc1_day(:,1) - taskAcc1_day(:,2);
plot(temp);plot([1 size(taskAcc1_day,1)],[0 0],'Color',[0.8 0.8 0.8]); title(num2str(nanmean(temp)))
subplot(1,3,2); hold on; temp =  taskAcc2_day(:,1) - taskAcc2_day(:,2);
plot(temp);plot([1 size(taskAcc1_day,1)],[0 0],'Color',[0.8 0.8 0.8]); title(num2str(nanmean(temp)))
subplot(1,3,3); hold on; temp = taskAcc2_day(:,1)+taskAcc1_day(:,1) - taskAcc1_day(:,2)-taskAcc2_day(:,2); 
plot(temp);plot([1 size(taskAcc1_day,1)],[0 0],'Color',[0.8 0.8 0.8]); title(num2str(nanmean(temp)))
%% PLOT 99 -- OLD WAY OF PLOTTING SIMPLE SCHEMATIC
acc_task1 = {}; acc_task2 = {}; acc_inter = {};
RT_task1 = {}; RT_task2 = {}; RT_inter = {};
trial_task1 = {2462:2501,2582:2621,2738:2777};
trial_task2 = {2343:2382, 2423:2461,2622:2659,2699:2737};
trial_inter = {2383:2422,2502:2581,2660:2698,2778:2855};
start = 2343; last = 2855;
ctxtFlag = zeros(last-start+1,1);
stimulus = zeros(last-start+1,1);

acc_combined = smoothdata(-(mouseObj.behav.responseType(start:last))+2,'movmean',30);

for i = 1:3
    acc_task1{i}  = -(mouseObj.behav.responseType(trial_task1{i}))+2;
    RT_task1{i}  = (mouseObj.behav.reactionTime(trial_task1{i}));
    ctxtFlag(trial_task1{i}-(start-1)) = 1;
    stimulus(trial_task1{i}-(start-1)) = mouseObj.behav.stimulus(trial_task1{i});
end

for i = 1:4
    acc_task2{i}  = -(mouseObj.behav.responseType(trial_task2{i}))+2;
    RT_task2{i}  = (mouseObj.behav.reactionTime(trial_task2{i}));
    ctxtFlag(trial_task2{i}-(start-1)) = 2;
    stimulus(trial_task2{i}-(start-1)) = mouseObj.behav.stimulus(trial_task2{i});
end

for i = 1:4
    acc_inter{i}  = -(mouseObj.behav.responseType(trial_inter{i}))+2;
    RT_inter{i}  = (mouseObj.behav.reactionTime(trial_inter{i}));
    ctxtFlag(trial_inter{i}-(start-1)) = 3;
    stimulus(trial_inter{i}-(start-1)) = mouseObj.behav.stimulus(trial_inter{i});
end

acc = [  mean(acc_task2{1}) mean(acc_inter{1}) ...
 mean(acc_task2{2})  mean(acc_task1{1}) mean(acc_inter{2})...
 mean(acc_task1{2})  mean(acc_task2{3}) mean(acc_inter{3})...
 mean(acc_task2{4})  mean(acc_task1{3}) mean(acc_inter{4}) ];


RT = [  mean(RT_task2{1}) mean(RT_inter{1}) ...
 mean(RT_task2{2})  mean(RT_task1{1}) mean(RT_inter{2})...
 mean(RT_task1{2})  mean(RT_task2{3}) mean(RT_inter{3})...
 mean(RT_task2{4})  mean(RT_task1{3}) mean(RT_inter{4}) ];

figure; plot(acc)

figure; hold on;
temp = acc_combined; temp(ctxtFlag~=1) = nan;
plot(temp,'Color',matlabColors(1),'LineWidth',3);
temp = acc_combined; temp(ctxtFlag~=2) = nan;
plot(temp,'Color',matlabColors(2),'LineWidth',3);
temp = acc_combined; temp(ctxtFlag~=3) = nan;
plot(temp,'Color',matlabColors(4),'LineWidth',3);
plot([0 length(temp)],[0.5 0.5],'Color',[0.8 0.8 0.8],'LineWidth',2);
xlim([0 length(temp)]); ylim([0.4 1]);

stimulus(stimulus==3 | stimulus==4) = stimulus(stimulus==3 | stimulus==4)+2;
figure; hold on;
temp = stimulus; temp(ctxtFlag~=1) = nan;
plot(temp,'.','Color',matlabColors(1),'LineWidth',2);
temp = stimulus; temp(ctxtFlag~=2) = nan;
plot(temp,'.','Color',matlabColors(2),'LineWidth',2);
temp = stimulus; temp(ctxtFlag~=3) = nan;
plot(temp,'.','Color',matlabColors(4),'LineWidth',2);
xlim([0 length(temp)]);

%% PLOT 3

clear; 
mouse ='zz101';
task = 'PT_pip_task2';

rootPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\';
loadPath = [rootPath filesep 'trialData\'];
load([loadPath mouse '.mat']);

mouseObj = wheel2AFC(trialData,'opsPath',['C:\Users\zzhu34\Documents\tempdata\octoData\opsFiles\' task ],'mouse',mouse );
mouseObj = mouseObj.removeMiss();  mouseObj = mouseObj.getAcc();

days = unique(mouseObj.behav.day); dayLen = [];

for i = 1:length(days); dayLen(i) = sum(mouseObj.behav.day == days(i)); end
temp = [0 cumsum(dayLen)]; temp = round((temp(1:end-1) + temp(2:end))/2);
load('C:\Users\zzhu34\Documents\tempdata\octoData\temp_task1_retention_data.mat','perf');
perf = perf{3}; 

perf = perf(:,1)./(perf(:,1)+perf(:,2))/2 + perf(:,4)./(perf(:,4)+perf(:,5))/2;

figure;
plot(mouseObj.behav.acc); hold on; plot(temp,perf,'--o')
%% FUNCTIONS 
function [behav, task1Idx,task2Idx,dayTransIdx] = selectInterLeaveDay(mouseObj,taskDay)
    interleaveDay = length(taskDay);
    
    interleaveIdx = (mouseObj.behav.day == taskDay(1));
    for i = 1:interleaveDay
        dayTransIdx(i) =  find(mouseObj.behav.day == taskDay(i),1,'last') - ...
            find(mouseObj.behav.day == taskDay(i),1,'first') + 1;
        interleaveIdx = interleaveIdx | mouseObj.behav.day == taskDay(i);
    end
    
    
    behav = mouseObj.behav(interleaveIdx,:);
    
    
    dayTransIdx = cumsum(dayTransIdx);
    task1Idx = (behav.stimulus==1 | behav.stimulus==2);
    task2Idx = (behav.stimulus==3 | behav.stimulus==4);
    
end


function [ctxtFlag, ctxtIdx,blockBin,blockCtxt, blockDay ]=getContextIdx(daylen)

dayLabel = {'interleave','20block','40block','randomblock'};
labelIdx = {1,6:17,[2:5 18:length(daylen)],[]};

daylen = [0 daylen];
ctxtFlag = ones(1,daylen(end));

for i = 1:length(dayLabel)

    for j = 1:length(labelIdx{i})
        tempIdx = labelIdx{i}(j);
        dayIdx = (daylen(tempIdx)+1):daylen(tempIdx+1);
        
        switch dayLabel{i}

        case 'interleave'
            tempCtxt = ones(1,length(dayIdx)) * 3;             
        case '40block'
            tempCtxt = [ones(1,40) ones(1,40)*2 ones(1,40)*3 ones(1,40)*2 ...
             ones(1,40) ones(1,80)*3];
        case '20block'
            tempCtxt = [ones(1,20) ones(1,20)*2 ones(1,20)*3 ones(1,20)*2 ...
             ones(1,20) ones(1,40)*3];
            tempCtxt = [tempCtxt tempCtxt];
        case 'randomblock'
            tempCtxt = ones(1,length(dayIdx)) * 4; 
        end 
        ctxtFlag(dayIdx) = tempCtxt;
    end 
end
blockBin = [0 find(diff(ctxtFlag)~=0)]; 
blockCtxt = zeros(1,length(blockBin)-1);
ctxtIdx = cell(3,1);blockDay = [];
for i = 1:length(blockBin)-1
    tempCtxt = ctxtFlag(blockBin(i)+1);
    if any(tempCtxt == [1 2 3])
        ctxtIdx{tempCtxt} = cat(1,ctxtIdx{tempCtxt},[blockBin(i)+1 blockBin(i+1)]);  
    end
    blockCtxt(i) = tempCtxt;
    
    blockDay(i) = find(blockBin(i)+1>daylen,1,'last');
end 
end

