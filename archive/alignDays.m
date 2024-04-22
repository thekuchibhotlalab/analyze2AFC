function trialData = alignDays(trialData)


for i = 1:length(trialData.stimulus)
    if ~isempty(trialData.wheelSoundOn{i})
        trialData.wheelSoundOn{i} = nansum(~isnan(trialData.wheelSoundOn{i}),2);
    else
        trialData.wheelSoundOn{i}=nan(length(trialData.stimulus{i}),1);
    end
    
    if ~trialData.wheelSoundOnCheckFlag{i}
        trialData.wheelPreSound{i} = nan(length(trialData.stimulus{i}),1);
        trialData.wheelSoundOn{i} = nan(length(trialData.stimulus{i}),1);
    end
    
end
trialDataKeepMiss = trialData;
rmFieldnames = {'date','trainingType','wheelSoundOnCheckFlag','sortIdx'};
for i = 1:length(rmFieldnames)
    if isfield(trialData,rmFieldnames{i})
        trialData = rmfield(trialData,rmFieldnames{i});
    end
end
learningCurveBin = 50;

% calculate trial each day 
dayLen = cellfun(@length, trialDataKeepMiss.action);
dayLenNoMiss = cellfun(@(x)(sum(x~=0)), trialDataKeepMiss.action);
lowActionDayFlag = zeros(sum(dayLenNoMiss),1);
tempActionRate = dayLenNoMiss./dayLen;
dayLenNoMissCumSum = [0 cumsum(dayLenNoMiss)];
for i = 1:length(dayLenNoMiss)
    if tempActionRate(i)<0.6
        lowActionDayFlag(dayLenNoMissCumSum(i)+1:dayLenNoMissCumSum(i+1)) = 1;
    end
end

% concatenate trials
trialData = structfun(@(x)fn_cell2mat(x,1),trialData,'UniformOutput',false);

missFlag = trialData.action==0;
probeFlag = trialData.context==3;

trialData = structfun(@(x)fn_removeIdx(x,missFlag),trialData,'UniformOutput',false);

trialData.accuracy = smoothdata(double(trialData.responseType==1),'movmean',learningCurveBin);

[trialData.bias, trialData.acc_L, trialData.acc_R] = fn_getBias(trialData.stimulus,trialData.responseType);

[probeData, trialNum,probeTrialNum,trialNumNoMiss,probeTrialNumNoMiss,reinfDataBef, reinfDataAft]...
    = cellfun(@fn_getProbe,trialDataKeepMiss.stimulus,trialDataKeepMiss.responseType,trialDataKeepMiss.context,'UniformOutput',false);

actL = trialData.action == 1;biasBin = 20;
actL = smoothdata(actL,'movmean',biasBin);
actAxis = biasBin/2+1:length(actL)-biasBin/2;

biasThreshold = 0.15; 
biasL = find(actL >= (biasThreshold+0.5) ); biasR = find(actL <= (0.5-biasThreshold));
trialBlockThreshold = 5;

biasL_incre = diff(biasL); temp = biasL_incre(2:end);
startFlag = find(biasL_incre>=5); 
trialData.biasBlockL_start = biasL([1; startFlag+1]) ; 
trialData.biasBlockL_end = biasL([startFlag; length(biasL)]) ;
trialData.biasBlockL_len= trialData.biasBlockL_end - trialData.biasBlockL_start+1;
trialData.biasBlockL_start = trialData.biasBlockL_start(trialData.biasBlockL_len>=trialBlockThreshold);
trialData.biasBlockL_end = trialData.biasBlockL_end(trialData.biasBlockL_len>=trialBlockThreshold);
trialData.biasBlockL_len = trialData.biasBlockL_len(trialData.biasBlockL_len>=trialBlockThreshold);

biasR_incre = diff(biasR); temp = biasR_incre(2:end);
startFlag = find(biasR_incre>=5); 
trialData.biasBlockR_start = biasR([1; startFlag+1]) ; 
trialData.biasBlockR_end = biasR([startFlag; length(biasR)]) ;
trialData.biasBlockR_len= trialData.biasBlockR_end - trialData.biasBlockR_start+1;
trialData.biasBlockR_start = trialData.biasBlockR_start(trialData.biasBlockR_len>=trialBlockThreshold);
trialData.biasBlockR_end = trialData.biasBlockR_end(trialData.biasBlockR_len>=trialBlockThreshold);
trialData.biasBlockR_len = trialData.biasBlockR_len(trialData.biasBlockR_len>=trialBlockThreshold);

trialData.probeData = fn_cell2mat(probeData,1);
trialData.reinfDataBef = fn_cell2mat(reinfDataBef,1);
trialData.reinfDataAft = fn_cell2mat(reinfDataAft,1);

trialNum = fn_cell2mat(trialNumNoMiss,1);
probeTrialNum = fn_cell2mat(probeTrialNumNoMiss,1);
trialData.probeDayNum = find(~isnan(probeTrialNum));

cumsumTrials = [0;cumsum(trialNum)];
probeTrialNum = probeTrialNum + cumsumTrials(1:end-1);
probeTrialNum(isnan(probeTrialNum)) = [];

trialData.trialNum = trialNum;
trialData.probeTrialNum = probeTrialNum;


trialData.reactionTime = trialData.reactionTime;

trialData.lowActionDayFlag = lowActionDayFlag;
end