classdef wheel2AFC < twoChoiceTask
    % wheel2AFC - sub-class of twoChoiceTask for ONE mouse
    %   wheel2AFC should contain all the training days from wheel training
    %   to 2AFC training
    %   
    
    properties
        % structure containing info about wheel movement
        probe
        biasBlock
    end
    
    methods
        %-----------------------------CONSTRUCTOR METHODS-----------------------------------------
        function obj = wheel2AFC(inStruct,varargin)
            %MOUSEOBJ Construct an instance of this class
            %   Detailed explanation goes here        
            p = fn_inputParser(); p.parse(varargin{:});
            inputNames = fieldnames(inStruct);
            
            % ---------- 1.1 Write variables into object ---------
            obj = obj@twoChoiceTask(inStruct,varargin{:},'subclassPause',p.Results.subclassPause);
            %obj = declareEmptyStruct(obj);
            % ---------- 1.2 Write input into MULTIDAY and WHEEL ---------
            %multiDayNames = {'date','trainingType','sortIdx'};
            %obj.multiday = write2struct(obj.multiday,multiDayNames);

            %wheelNames = {'wheelSoundOn','wheelPreSound','wheelSoundOnCheckFlag'};
            %obj.wheel = write2struct(obj.wheel,wheelNames);           
            % ---------- 2.0  -----------------------
            %if p.Results.directImport
            %    taskObjPropName = properties('wheel2AFC');
            %    % WRITE variables into OBJ
            %    for i = 1:length(taskObjPropName)
            %        if contains(taskObjPropName{i},inputNames) ;obj.(taskObjPropName{i}) = inStruct.(taskObjPropName{i}); end
            %    end
            %    return;
            %end      
            % ---------- 3.0 Select Protocol and extract probe if is cell ---------
            %if strcmp(obj.ops.inputType,'cell')
            %    if ~isempty(p.Results.selectProtocol); obj = selectProtocol(obj,p.Results.selectProtocol); end
            %    % Correct the number of wheel trials
            %    if ~isempty(p.Results.correctWheelTrial); obj = correctWheelTrial(obj); end
            %end
            % ---------- 3.1 Concatenate across days --------   
            %if strcmp(obj.ops.inputType,'cell') && strcmp(obj.ops.multidayType,'mat')
            %    % Construct probe flag with information of each day
            %    [obj.probeFlag.probeFlag, obj.probeFlag.probeStartFlag, obj.probeFlag.reinfBefFlag, obj.probeFlag.reinfAftFlag] = cellfun(@fn_getProbeFlag,obj.context,num2cell(1:length(obj.stimulus)),'UniformOutput',false);
            %    % Concatenate trials across days
            %    obj = concatenateDay(obj);
            %    obj = separateProbe (obj);  
            %end
            % ---------- 3.2 Remove miss trials, if instructed-----------------------
            %if obj.ops.removeMiss; missFlag = (obj.action == 0); obj = removeIdx(obj,missFlag); end

            %obj = getAcc(obj);
            %obj = getChoiceAcc(obj);
            %obj = getBiasBlock(obj);
            
            % wheel trajectory code
            if contains('wheelPos_aligned',fields(inStruct)) && size(inStruct.wheelPos_aligned{1},1) == 1
                for j = 1:length(inStruct.wheelPos_aligned); inStruct.wheelPos_aligned{j} = inStruct.wheelPos_aligned{j}'; end
            end

            if ~obj.ops.preprocessedInput; inStruct = fn_preprocessInput(inStruct,obj.ops); end
            if contains('wheelPos_aligned',fields(inStruct)) && ~isempty(inStruct.wheelPos_aligned); processWheel(inStruct); 
            else; disp([obj.ops.mouse ' -- no Wheel position data!']); end
            %catch; disp([obj.ops.mouse ' -- no Wheel Timestamp Data'])
            %end
            
            function processWheel(inStruct)
                timeBef = 3000; timeAft = 5000; % specify the time point that is taken for plotting
                temp = inStruct.wheelPos_aligned; 
                tempCrop = cellfun(@correctEmptyWheelTrace,temp,'UniformOutput',false);
                obj.behavVar.wheelPos_aligned = fn_cell2matFast(tempCrop,1);
                %obj.behavVar.wheelPos_aligned = inStruct.wheelPos_aligned;
                outStruct = fn_preprocessWheel(temp,obj.behavVar.wheelPos_aligned,timeBef,timeAft,obj.behav,obj.ops);
                a = fields(outStruct); for k = 1:length(a); obj.behav.(a{k}) = outStruct.(a{k}); end
                function x = correctEmptyWheelTrace(x)
                    if isempty(x); x = zeros(1,timeBef+timeAft+1);
                    else; x = x(1:timeBef+timeAft+1);
                    end

                end

            end

            

            function p = fn_inputParser()
                p = inputParser;
                p.KeepUnmatched = true;
                % ------ Varargin for superclass
                p.addParameter('opsPath', pwd);
                p.addParameter('mouse', []);
                p.addParameter('verbose', true);
                p.addParameter('subclassPause', true);
                % ------ Varargin for this class

                p.addParameter('directImport', false);
            end
            
            function S = write2struct(S,varNames)
                for i = 1:length(varNames)
                    if contains(varNames{i},inputNames)
                        S.(varNames{i}) = inStruct.(varNames{i});
                    end 
                end  
            end
            
            %function obj = declareEmptyStruct(obj)
            %    obj.probe = struct(); obj.probeFlag = struct();
            %end
            
            
        end
        %-----------------------------METHODS EXTENDING FROM TWOCHOICETASK CLASS-----------------------------------------
%         function obj = selectDay(obj,selectFlag)
%             obj = selectDay@twoChoiceTask(obj,selectFlag);
%             multidayPropSel = {'date','trainingType','sortIdx'};
%             for i = 1:length(multidayPropSel); obj.multiday.(multidayPropSel{i}) = obj.multiday.(multidayPropSel{i})(selectFlag);end
%             wheelPropSel = {'wheelPreSound','wheelSoundOn','wheelSoundOnCheckFlag'};
%             for i = 1:length(wheelPropSel); obj.wheel.(wheelPropSel{i}) = obj.wheel.(wheelPropSel{i})(selectFlag);end
%         end

%         function obj = removeIdx(obj,idx)
%             obj = removeIdx@twoChoiceTask(obj,idx);
%             obj.wheel = fn_structfun(@fn_removeIdx, obj.wheel,'opsArgIn',{idx},'funcName', 'REMOVE_IDX','verbose', obj.ops.verbose);
%             obj.probe = fn_structfun(@fn_removeIdx, obj.probe,'opsArgIn',{idx},'funcName', 'REMOVE_IDX','verbose', obj.ops.verbose);
%             obj.probeFlag = fn_structfun(@fn_removeIdx, obj.probeFlag,'opsArgIn',{idx},'funcName', 'REMOVE_IDX','verbose', obj.ops.verbose);
%         end
        
        
        %-----------------------------METHODS SPECIFIC TO WHEEL 2AFC CLASS-----------------------------------------
        function obj = calculateProbe(obj) % this is the old method of calculating and binning probe
            catDay = false;
            probeDays = unique(obj.behav.day(logical(obj.behav.goodProbe)));
            
            obj.probe.probeData = zeros(length(probeDays),6);
            obj.probe.befData = zeros(length(probeDays),6);
            obj.probe.aftData = zeros(length(probeDays),6);
            obj.probe.trialIdx = zeros(length(probeDays),1);
            
            if catDay
                for i = 1:length(probeDays)
                    probeIdx = (obj.behav.day==probeDays(i)) & (obj.behav.goodProbe==1);
                    befIdx = (obj.behav.day==probeDays(i)) & (obj.behav.reinfBef==1);
                    aftIdx = (obj.behav.day==probeDays(i)) & (obj.behav.reinfAft==1);

                    [probeOn, probeOff] = fn_getBlockOnOff(probeIdx);

                    obj.probe.probeData(i,:) = fn_getIdxAccBias(obj.behav.stimulus(probeIdx), obj.behav.responseType(probeIdx));
                    obj.probe.befData(i,:) = fn_getIdxAccBias(obj.behav.stimulus(befIdx), obj.behav.responseType(befIdx));
                    obj.probe.aftData(i,:) = fn_getIdxAccBias(obj.behav.stimulus(aftIdx), obj.behav.responseType(aftIdx));
                    obj.probe.trialIdx(i) = mean(find(probeOn));
                    disp([sum(probeIdx) sum(befIdx) sum(aftIdx)])

                end
            else   
                % REDO the bef and aft selection since miss is removed
                [obj.probe.onIdx, obj.probe.offIdx, probeIdx] = fn_getBlockOnOff(obj.behav.goodProbe==1);
                obj.behav.reinfBef = zeros(size(obj.behav.reinfBef)); obj.behav.reinfAft = zeros(size(obj.behav.reinfAft));
                obj.probe.day = [];
                befAftLen = 30;
                for i = 1:length(probeIdx)
                    tempProbeIdx = probeIdx{i}; tempDay = obj.behav.day(tempProbeIdx(1)); tempDayEnd = find(obj.behav.day==tempDay,1,'last');
                    obj.probe.probeData(i,:) = fn_getTrialTypes(obj.behav.stimulus(tempProbeIdx), obj.behav.responseType(tempProbeIdx));
                    % old way of taking 10 trials pre and post
                    %tempBefIdx = tempProbeIdx - length(tempProbeIdx); tempAftIdx = tempProbeIdx + length(tempProbeIdx);
                    % new way of calculating pre and post
                    tempBefIdx = tempProbeIdx(1)-befAftLen-1:tempProbeIdx(1)-2; % use -2 instead of -1 to take out the 1 correct trial before
                    tempAftIdx = tempProbeIdx(end)+1:tempProbeIdx(end)+befAftLen; 
                    tempFirstRewardIdx = tempAftIdx(find(obj.behav.reward(tempAftIdx)==1,1)); excludeLen = 10;
                    tempAftIdx = tempFirstRewardIdx+1+excludeLen:tempFirstRewardIdx+befAftLen+excludeLen;
                    tempAftIdx(tempAftIdx>tempDayEnd) = [];

                    obj.behav.reinfBef(tempBefIdx) = 1; obj.behav.reinfAft(tempAftIdx) = 1; 
                    obj.behav.badReinf(tempProbeIdx(1)-1) = 1; 
                    if ~isempty(tempFirstRewardIdx); obj.behav.badReinf(tempProbeIdx(end)+1:tempFirstRewardIdx) = 1; else; disp(['Empty idx after probe, day=' int2str(tempDay)]); end
                    obj.probe.befData(i,:) = fn_getTrialTypes(obj.behav.stimulus(tempBefIdx), obj.behav.responseType(tempBefIdx));
                    obj.probe.aftData(i,:) = fn_getTrialTypes(obj.behav.stimulus(tempAftIdx), obj.behav.responseType(tempAftIdx));
                    obj.probe.befAftData(i,:) = fn_getTrialTypes(obj.behav.stimulus([tempBefIdx tempAftIdx]), obj.behav.responseType([tempBefIdx tempAftIdx]));

                    nProbe = length(tempProbeIdx); tempFirstIdx = tempProbeIdx(1:floor(nProbe/2)); tempLastIdx = tempProbeIdx(floor(nProbe/2)+1:end); 
                    obj.probe.firstHalfData(i,:) = fn_getTrialTypes(obj.behav.stimulus(tempFirstIdx), obj.behav.responseType(tempFirstIdx));
                    obj.probe.lastHalfData(i,:) = fn_getTrialTypes(obj.behav.stimulus(tempLastIdx), obj.behav.responseType(tempLastIdx));
                    
                    catchResp = obj.behav.action(obj.probe.onIdx(i)-1); if catchResp==2; catchResp=-1; end 
                    obj.probe.catch(i) = catchResp; [tempBias] = fn_getAccBiasByCount(obj.probe.befData(i,:));
                    obj.probe.catchBias(i) = catchResp*tempBias;
                    obj.probe.idx = probeIdx;
                    obj.probe.day(i) = tempDay;
                    if isfieldTable( obj.behav, 'modelBias')
                        obj.probe.befModelBias(i) = mean(obj.behav.modelBias(tempBefIdx));
                        obj.probe.aftModelBias(i) = mean(obj.behav.modelBias(tempAftIdx));
                    else
                        obj.probe.befModelBias(i) = nan;
                        obj.probe.aftModelBias(i) = nan;
                    end
                end
            end
        end
        
        % new methods for determining the probe bin
        function [probeBin] = binProbe(obj,binFlag) % combine adjacent bins to produce enough count
            if nargin == 1; binFlag  = true; end 
            nProbeBlock  = size(obj.probe.probeData,1);
            if binFlag
                probeBin = struct();
                trialThre = 8;
                tempProbeBin = zeros(1,6); binIdx = [];
                for i = 1:nProbeBlock
                    tempProbeBin = tempProbeBin + obj.probe.probeData(i,:);
                    if (sum(tempProbeBin(1:2)) > trialThre && sum(tempProbeBin(4:5)) > trialThre)
                        binIdx = [binIdx i];tempProbeBin = zeros(1,6);
                    end
                end
            else; binIdx = 1:nProbeBlock; end
            binIdx = [0 binIdx]; 
            if binIdx(end) > size(obj.probe.probeData,1); binIdx(end) = nProbeBlock; end 
            for i = 1:length(binIdx)-1
                tempIdx = (binIdx(i)+1):binIdx(i+1);
                [probeBin.probeBias(i,1),probeBin.probeAcc(i,1)] = fn_getAccBiasByCount(sum(obj.probe.probeData(tempIdx,:),1));
                [probeBin.befBias(i,1),probeBin.befAcc(i,1)] = fn_getAccBiasByCount(sum(obj.probe.befData(tempIdx,:),1));
                [probeBin.aftBias(i,1),probeBin.aftAcc(i,1)] = fn_getAccBiasByCount(sum(obj.probe.aftData(tempIdx,:),1));
                [probeBin.firstHalfBias(i,1),probeBin.firstHalfAcc(i,1)] = fn_getAccBiasByCount(sum(obj.probe.firstHalfData(tempIdx,:),1));
                [probeBin.lastHalfBias(i,1),probeBin.lastHalfAcc(i,1)] = fn_getAccBiasByCount(sum(obj.probe.lastHalfData(tempIdx,:),1));
                

                probeBin.befModelBias(i,1) = nanmean(obj.probe.aftModelBias(tempIdx));
                probeBin.aftModelBias(i,1) = nanmean(obj.probe.aftModelBias(tempIdx));
                probeBin.onIdx(i,1) = obj.probe.onIdx(tempIdx(1)); probeBin.offIdx(i,1) = obj.probe.offIdx(tempIdx(end));
                probeBin.probeIdx(i,1) = (probeBin.onIdx(i,1)+ probeBin.offIdx(i,1))/2;
            end
            
        end 
        function [probeBin] = binProbeByDay(obj)
            uniqueDay = unique(obj.probe.day);
            for i = 1:length(uniqueDay)
                tempIdx = (uniqueDay(i) == obj.probe.day);
                [probeBin.bias(i),probeBin.acc(i)] = fn_getAccBiasByCount(sum(obj.probe.probeData(tempIdx,:),1));
                probeBin.trial(i) = round(mean(cell2mat(obj.probe.idx(tempIdx))));
                probeBin.nTrial(i) = length(cell2mat(obj.probe.idx(tempIdx)));
                probeBin.day(i) = uniqueDay(i); 
            end
        end
        % VISUALIZATION -- BIN PROBE TRIAL ACCORDING TO GIVEN TRIAL WINDOW 
        function [probeBin] = binProbeByTrial(obj,trialBinStart,trialBinEnd)
            nBins = length(trialBinStart) ;
            %[~,~,binIdx] = histcounts(obj.probe.onIdx,trialBin);
            binIdx = zeros(length(obj.probe.onIdx),nBins);
            for i = 1:length(obj.probe.onIdx)
                tempIdx = (obj.probe.onIdx(i) >= trialBinStart) & (obj.probe.onIdx(i) < trialBinEnd);
                binIdx(i,:) = tempIdx;
            end

            for i = 1:nBins
                tempIdx = logical(binIdx (:,i));
                probeBin.probeData(i,:) = nansum(obj.probe.probeData(tempIdx,:),1);
                probeBin.befData(i,:) = nansum(obj.probe.befData(tempIdx,:),1);
                probeBin.aftData(i,:) = nansum(obj.probe.aftData(tempIdx,:),1);
                probeBin.befAftData(i,:) = nansum(obj.probe.befAftData(tempIdx,:),1);

                [probeBin.probeBias(i,1),probeBin.probeAcc(i,1),probeBin.probeDp(i,1)] = fn_getAccBiasByCount(nansum(obj.probe.probeData(tempIdx,:),1));
                [probeBin.befBias(i,1),probeBin.befAcc(i,1),probeBin.befDp(i,1)] = fn_getAccBiasByCount(nansum(obj.probe.befData(tempIdx,:),1));
                [probeBin.aftBias(i,1),probeBin.aftAcc(i,1),probeBin.aftDp(i,1)] = fn_getAccBiasByCount(nansum(obj.probe.aftData(tempIdx,:),1));
                [probeBin.firstHalfBias(i,1),probeBin.firstHalfAcc(i,1),probeBin.firstHalfDp(i,1)] = fn_getAccBiasByCount(nansum(obj.probe.firstHalfData(tempIdx,:),1));
                [probeBin.lastHalfBias(i,1),probeBin.lastHalfAcc(i,1),probeBin.lastHalfDp(i,1)] = fn_getAccBiasByCount(nansum(obj.probe.lastHalfData(tempIdx,:),1));
                
                
                probeBin.befModelBias(i,1) = nanmean(obj.probe.aftModelBias(tempIdx));
                probeBin.aftModelBias(i,1) = nanmean(obj.probe.aftModelBias(tempIdx));
            end   
        end
        function [probeData, trialBin,learningOnset, onsetBin] = binProbeByTrialFromLearningOnset(obj,accThre,windowSize)
            [trialBin, learningOnset] = obj.fn_findLearningOnset(accThre,windowSize);
            windowStart = trialBin(1:end-1); windowEnd = trialBin(2:end);
            probeData = obj.binProbeByTrial(windowStart,windowEnd);
            onsetBin = find(learningOnset < trialBin, 1, 'first')-1;
        end
        % REINF VS. PROBE COMPARISON, (BEF VS. AFT) OR (SUBSAMPLE WITHIN DAY)
        function [probeDataAllTrial,probeByDay,nProbePerDay,probeByDayBin,trialLim,probeBlockID] = computeProbeAlignByTrial(obj,compType,accThre,windowSize)
            if nargin == 1; accThre = [0.7 nan]; windowSize= 400; end
            [~, learningOnset] = obj.fn_findLearningOnset(accThre,windowSize);
            selTrial = learningOnset-700:learningOnset+700;
            nTrials = size(obj.behav,1); 
            selTrial(selTrial<=0) = []; selTrial(selTrial>nTrials) = [];
            dayStart = obj.behav.day(selTrial(1)); dayEnd = obj.behav.day(selTrial(end));
            disp([obj.ops.mouse ' start day ' int2str(dayStart) '; end day ' int2str(dayEnd) '; total day=' int2str(dayEnd-dayStart+1)])
            if strcmp(compType,'befAft')
                [probeByDay,nProbePerDay,probeByDayBin,trialLim] = fn_combineProbeByDay(obj,[dayStart dayEnd],compType);
                [probeDataAllTrial] = fn_combineProbeAllTrial(obj,[dayStart dayEnd],compType);
                [probeBlockID] = fn_probeBlockID(obj,[dayStart dayEnd]);
            elseif strcmp(compType,'day')
                [probeByDay,nProbePerDay,probeByDayBin,trialLim] = fn_combineProbeByDay(obj,[dayStart dayEnd],compType);
                [probeDataAllTrial] = fn_combineProbeAllTrial(obj,[dayStart dayEnd],compType);
                probeBlockID = struct(); % do not evaluate 
            end
            
            nDay = length(dayStart:dayEnd);
            trialLim(1) = find(obj.behav.day==dayStart,1); trialLim(2) = find(obj.behav.day==dayEnd,1,'last'); %nBefAft = 30;


        end
        % REINF VS. PROBE COMPARISON, randomized control
        function [probeRand,probeRand2] = computeProbeAlignByTrialRandomize(obj)
            accThre = [0.7 nan]; windowSize= 400; 
            [~, learningOnset] = obj.fn_findLearningOnset(accThre,windowSize);
            selTrial = learningOnset-300:learningOnset+300;
            nTrials = size(obj.behav,1); 
            selTrial(selTrial<=0) = []; selTrial(selTrial>nTrials) = [];
            dayStart = obj.behav.day(selTrial(1)); dayEnd = obj.behav.day(selTrial(end));
            nDay = length(dayStart:dayEnd);
            trialLim(1) = find(obj.behav.day==dayStart,1); trialLim(2) = find(obj.behav.day==dayEnd,1,'last'); %nBefAft = 30;
            % COMPUTE BIAS BY EACH DAY (by subsampling trials)
            probeRand = cell(1,nDay); probeRand2 = cell(1,nDay); 
            for i = 1:nDay
                tempDay = dayStart+i-1; 
                tempBehav = obj.behav(obj.behav.day == tempDay,:); 
                reinfFlag = ~tempBehav.goodProbe & ~tempBehav.badReinf; 
                reinfIdx = find(reinfFlag);

                nProbe = sum(tempBehav.goodProbe);
                tempBin = 0:floor(nProbe/2):(length(reinfIdx)-nProbe);
                [tempBiasBase,tempAccBase,~,~] = fn_subsample(tempBehav,find(reinfFlag),1000,nProbe); 
                tempBiasBase = nanmean(abs(tempBiasBase)); tempAccBase = nanmean(abs(tempAccBase));
                sampleRandomized = [];sampleRandomized2 = [];
                for j = 1:length(tempBin)
                    tempIdx = reinfIdx(tempBin(j)+1:tempBin(j)+nProbe);
                    [tempBias, tempAcc, tempDP, tempCri, ~] = fn_getAccBias(tempBehav.stimulus(tempIdx), ...
                        tempBehav.responseType(tempIdx)==1,tempBehav.responseType(tempIdx)==0);
                    sampleRandomized(:,j) = [abs(tempBiasBase) - abs(tempBias) tempAcc - tempAccBase];
                    sampleRandomized2(:,j) = [abs(tempBiasBase) abs(tempBiasBase)-abs(tempBias) ];
                end

                
                %nProbe = sum(~reinfFlag); nRep = 100;
                %sampleRandomized = zeros(2,nRep); sampleRandomized2 = zeros(2,nRep);
                %for j = 1: nRep
                    %[tempBias1,tempAcc1,~,~] = fn_subsample(tempBehav,find(reinfFlag),1,nProbe); 
                    %[tempBias2,tempAcc2,~,~] = fn_subsample(tempBehav,find(reinfFlag),1,nProbe); 
                    %sampleRandomized(:,j) = [abs(tempBias2) - abs(tempBias1) tempAcc1 - tempAcc2];
                    %sampleRandomized2(:,j) = [abs(tempBias2) abs(tempBias2) - abs(tempBias1) ];
                %end
                probeRand{i} = sampleRandomized; probeRand2{i} = sampleRandomized2;
            end
            probeRand = fn_cell2mat(probeRand,2); probeRand2 = fn_cell2mat(probeRand2,2);
        end
        
        function [probeDataAllTrial,probeByDay,nProbePerDay,probeByDayBin,trialLim] = computeProbeAlignByPerf(obj,compType,selThre,reinfAccBin)
            if nargin == 1; selThre = [0.6 0.8]; reinfAccBin= 400; end
            nTrials = size(obj.behav,1); 
            [~,reinfAcc] = fn_getAccBiasSmooth(obj.behav.stimulus,obj.behav.responseType,reinfAccBin);
            reinfStartIdx = find(reinfAcc > selThre(1),1);
            reinfEndIdx = find(reinfAcc > selThre(2),1); reinfEndIdx = min([nTrials,reinfEndIdx]);

            dayStart = obj.behav.day(reinfStartIdx); dayEnd = obj.behav.day(reinfEndIdx);
            if strcmp(compType,'befAft')
                [probeByDay,nProbePerDay,probeByDayBin,trialLim] = fn_combineProbeByDay(obj,[dayStart dayEnd],compType);
                [probeDataAllTrial] = fn_combineProbeAllTrial(obj,[dayStart dayEnd]);
            elseif strcmp(compType,'day')
                [probeByDay,nProbePerDay,probeByDayBin,trialLim] = fn_combineProbeByDay(obj,[dayStart dayEnd],compType);
                [probeDataAllTrial] = fn_combineProbeAllTrial(obj,[dayStart dayEnd]);
            end
        end
        function [probeBlockID] = fn_probeBlockID(obj,selDay)
            dayStart = selDay(1); dayEnd = selDay(2);
            trialLim(1) = find(obj.behav.day==dayStart,1); trialLim(2) = find(obj.behav.day==dayEnd,1,'last'); 
            probeSelFlag = find(obj.probe.onIdx >= trialLim(1) & obj.probe.onIdx <= trialLim(2));
            reinfCount = zeros(1,6); probeCount = zeros(1,6);
            biasThre = 0.4;
            probeBlockID = [];
            for i = 1:length(probeSelFlag)
                tempIdx = probeSelFlag(i);
                tempProbe = obj.probe.probeData(tempIdx,:);
                tempBefBehav = obj.behav(obj.probe.onIdx(i)-31:obj.probe.onIdx(i)-2,:); %evaluate bias in 30 trials before probe onset
                [tempBias] = fn_getAccBias(tempBefBehav.stimulus,tempBefBehav.reward==1);
                catchResp = obj.behav.responseType(obj.probe.onIdx(i)-1);

                tempReinf = obj.probe.befAftData(tempIdx,:);
                tempProbe = obj.probe.probeData(tempIdx,:);

            end

        end

        function [probeData] = fn_combineProbeAllTrial(obj,selDay,compType)
            dayStart = selDay(1); dayEnd = selDay(2);

            trialLim(1) = find(obj.behav.day==dayStart,1); trialLim(2) = find(obj.behav.day==dayEnd,1,'last'); 
            probeSelFlag = find(obj.probe.onIdx >= trialLim(1) & obj.probe.onIdx <= trialLim(2));

            %probeSelFlag = find(obj.probe.onIdx >= reinfStartIdx & obj.probe.onIdx <= reinfEndIdx);

            probeCount = zeros(1,6);probeCountFirstHalf = zeros(1,6); probeCountLastHalf = zeros(1,6); 
            befCount = zeros(1,6); aftCount = zeros(1,6); 
            probeCountStrongViolation = zeros(1,6); probeCountWeakViolation = zeros(1,6); 
            for i = 1:length(probeSelFlag)
                tempIdx = probeSelFlag(i);
                tempProbe = obj.probe.probeData(tempIdx,:);
                tempProbeFirstHalf = obj.probe.firstHalfData(tempIdx,:);
                tempProbeLastHalf = obj.probe.lastHalfData(tempIdx,:);
                tempBef = obj.probe.befData(tempIdx,:);
                tempAft = obj.probe.aftData(tempIdx,:);
                %tempBefAft = obj.probe.befData(tempIdx,:) + obj.probe.aftData(tempIdx,:);
                %[tempBias,~,~] = fn_getAccBiasByCount(tempBefAft);
                %if tempBias<-0.2
                %    temp = tempProbe; temp(1:3) = tempProbe(4:6); temp(4:6) = tempProbe(1:3); tempProbe = temp; 
                %    temp = tempBef; temp(1:3) = tempBef(4:6); temp(4:6) = tempBef(1:3); tempBef = temp; 
                %    temp = tempAft; temp(1:3) = tempAft(4:6); temp(4:6) = tempAft(1:3); tempAft = temp; 
                %end
                probeCount = probeCount + tempProbe;
                probeCountFirstHalf = probeCountFirstHalf + tempProbeFirstHalf;
                probeCountLastHalf = probeCountLastHalf + tempProbeLastHalf;
                befCount = befCount + tempBef;
                aftCount = aftCount + tempAft;
                if obj.probe.catchBias(i)>0
                    probeCountStrongViolation = probeCountStrongViolation + tempProbe;
                else
                    probeCountWeakViolation = probeCountWeakViolation + tempProbe;
                end
            end

            nSample = 2000; 
            if strcmp(compType,'befAft')
                befAftCount = befCount + aftCount;
                
                [befBias,befAcc] = fn_calculateBiasBySample(befCount,nSample,[sum(probeCount(1:2)) sum(probeCount(3:4))]);
                [aftBias,aftAcc] = fn_calculateBiasBySample(aftCount,nSample,[sum(probeCount(1:2)) sum(probeCount(3:4))]);
                [befAftBias,befAftAcc] = fn_calculateBiasBySample(befAftCount,nSample,[sum(probeCount(1:2)) sum(probeCount(3:4))]);
                probeData.befBias = nanmean(abs(befBias)); probeData.befAcc = nanmean(abs(befAcc));
                probeData.aftBias = nanmean(abs(aftBias)); probeData.aftAcc = nanmean(abs(aftAcc));
                probeData.befAftBias = nanmean(abs(befAftBias)); probeData.befAftAcc = nanmean(abs(befAftAcc));
            elseif strcmp(compType,'day')
                tempBehav = obj.behav(trialLim(1):trialLim(2),:); 
                tempReinfIdx = find(~tempBehav.goodProbe & ~tempBehav.badReinf);
                [tempBias,tempAcc, tempDP, tempCri] = fn_subsample(tempBehav,tempReinfIdx,nSample,sum(probeCount));
                probeData.reinfBias = nanmean(abs(tempBias)); probeData.reinfAcc = nanmean(abs(tempAcc));
                probeData.reinfCri = nanmean(abs(tempCri)); probeData.reinfDP = nanmean(abs(tempDP));
            end 
            [probeData.probeBias,probeData.probeAcc,probeData.probeDP, probeData.probeCri] = fn_getAccBiasByCount(probeCount); 
            [probeData.probeBiasFirstHalf,probeData.probeAccFirstHalf,~] = fn_getAccBiasByCount(probeCountFirstHalf); 
            [probeData.probeBiasLastHalf,probeData.probeAccLastHalf,~] = fn_getAccBiasByCount(probeCountLastHalf); 
            [probeData.probeBiasStrong,probeData.probeAccStrong] = fn_getAccBiasByCount(probeCountStrongViolation); 
            [probeData.probeBiasWeak,probeData.probeAccWeak] = fn_getAccBiasByCount(probeCountWeakViolation); 
        end 
        function [probeByDay,nProbePerDay,probeByDayBin,trialLim] = fn_combineProbeByDay(obj,selDay,compType)
            dayStart = selDay(1); dayEnd = selDay(2);nDay = length(dayStart:dayEnd);
            trialLim(1) = find(obj.behav.day==dayStart,1); trialLim(2) = find(obj.behav.day==dayEnd,1,'last'); %nBefAft = 30;
            % COMPUTE BIAS BY EACH DAY (by subsampling trials)
            probeByDay.probe = nan(nDay,4); probeByDay.reinf = nan(nDay,4);
            probeByDay.bef = nan(nDay,4); probeByDay.aft = nan(nDay,4); 
            probeByDay.modelRemoveBias = nan(nDay,1); 
            for i = 1:nDay
                tempDay = dayStart+i-1; 
                tempBehav = obj.behav(obj.behav.day == tempDay,:); 
                probeFlag = tempBehav.goodProbe; 
                [probeByDay.probe(i,2), probeByDay.probe(i,1), probeByDay.probe(i,3), probeByDay.probe(i,4), ~] = fn_getAccBias(tempBehav.stimulus(probeFlag), ...
                    tempBehav.responseType(probeFlag)==1,tempBehav.responseType(probeFlag)==0);
                nProbe = sum(probeFlag); nRep = 5000;
                probeByDay.modelRemoveBias(i) = nanmean(tempBehav.modelPred_removeAllAcc(tempBehav.goodProbe));
                if strcmp(compType,'befAft') 
                    tempBefIdx = find(tempBehav.reinfBef); tempAftIdx = find(tempBehav.reinfAft);
                    if ~isempty(tempBefIdx) && length(tempBefIdx) >= nProbe                        
                        [tempBias,tempAcc,tempDP,tempCri] = fn_subsample(tempBehav,tempBefIdx,nRep,nProbe); 
                        probeByDay.bef(i,1) = nanmean(tempAcc); probeByDay.bef(i,2) = nanmean(abs(tempBias));
                        probeByDay.reinf(i,3) = nanmean(tempDP); probeByDay.reinf(i,4) = nanmean(abs(tempCri));
                    end 
                    if ~isempty(tempAftIdx) && length(tempAftIdx) >= nProbe
                        [tempBias,tempAcc,tempDP,tempCri] = fn_subsample(tempBehav,tempAftIdx,nRep,nProbe); 
                        probeByDay.aft(i,1) = nanmean(tempAcc); probeByDay.aft(i,2) = nanmean(abs(tempBias));
                        probeByDay.reinf(i,3) = nanmean(tempDP); probeByDay.reinf(i,4) = nanmean(abs(tempCri));
                    end 
                elseif strcmp(compType,'day')
                    tempReinfIdx = find(~tempBehav.badReinf & ~tempBehav.goodProbe);
                    [tempBias,tempAcc,tempDP,tempCri] = fn_subsample(tempBehav,tempReinfIdx,nRep,nProbe); 
                    %[tempBias,tempAcc] = fn_subsampleBlock(tempBehav,tempReinfIdx,nProbe); 
                    probeByDay.reinf(i,1) = nanmean(tempAcc); probeByDay.reinf(i,2) = nanmean(abs(tempBias));
                    probeByDay.reinf(i,3) = nanmean(tempDP); probeByDay.reinf(i,4) = nanmean(abs(tempCri));
                end
            end
            % BIN MULTIPLE DAYS TO MAKE SURE EACH DAY HAVE ENOUGH TRIALS (by subsampling trials)
            binThreshold = 14;
            nProbePerDay = [];
            for i = 1:nDay
                tempDay = dayStart+i-1; 
                tempBehav = obj.behav(obj.behav.day == tempDay,:); 
                probeFlag = tempBehav.goodProbe;
                nProbePerDay(i) = sum(probeFlag);
            end
            tempSum = 0; groupDay = {}; groupDay{1} = [];count = 1; 
            for i = 1:nDay
                tempDay = dayStart+i-1; 
                tempSum = tempSum+nProbePerDay(i); 
                groupDay{count} = [groupDay{count} tempDay];
                if tempSum>=binThreshold; tempSum = 0; count = count+1; groupDay{count} = [];end
            end
            if isempty(groupDay{end}); groupDay(end) = [];end
            if sum(nProbePerDay(groupDay{end}-dayStart+1))<=binThreshold; groupDay{end-1} = [groupDay{end-1} groupDay{end}]; groupDay(end) = [];end

            probeByDayBin.probe = nan(length(groupDay),2); probeByDayBin.bef = nan(length(groupDay),2); probeByDayBin.aft = nan(length(groupDay),2); 
            for i = 1:length(groupDay)
                tempDay = groupDay{i}; tempFlag = sum(obj.behav.day == tempDay,2)>0;
                tempBehav = obj.behav(tempFlag,:); 
                probeFlag = tempBehav.goodProbe; probeIdx = find(probeFlag);
                [probeByDayBin.probe(i,2), probeByDayBin.probe(i,1), ~, ~, ~] = fn_getAccBias(tempBehav.stimulus(probeFlag), ...
                            tempBehav.responseType(probeFlag)==1,tempBehav.responseType(probeFlag)==0);
                %{
                [onIdx, offIdx, ~, ~] = fn_getBlockOnOff(probeFlag);
                tempBefIdx = []; tempAftIdx = []; 
                for j = 1:length(onIdx)
                    
                    tempBefIdx = [tempBefIdx onIdx(j)-nBefAft:onIdx(j)-1]; tempAftIdx = [tempAftIdx offIdx(j)+1:offIdx(j)+nBefAft];
                    tempBefIdx(tempBefIdx<0) = []; tempBefIdx(tempBefIdx>length(probeFlag)) = []; 
                    tempFlag = sum(tempBefIdx==probeIdx,1); tempBefIdx(tempFlag>0) = [];
                    tempAftIdx(tempAftIdx<0) = []; tempAftIdx(tempAftIdx>length(probeFlag)) = []; 
                    tempFlag = sum(tempAftIdx==probeIdx,1); tempAftIdx(tempFlag>0) = []; 
                end
                %}
                % NEW method -- use the bef trials identified by previous codes
                tempBefIdx = find(tempBehav.reinfBef); tempAftIdx = find(tempBehav.reinfAft);

                nProbe = sum(probeFlag); nRep = 2000;
                if ~isempty(tempBefIdx) && length(tempBefIdx) >= nProbe
                    [tempBefBias,tempBefAcc] = fn_subsample(tempBehav,tempBefIdx,nRep,nProbe);
                    probeByDayBin.bef(i,1) = nanmean(tempBefAcc); probeByDayBin.bef(i,2) = nanmean(abs(tempBefBias));
                end 
                if ~isempty(tempBefIdx) && length(tempAftIdx) >= nProbe
                    [tempAftBias,tempAftAcc] = fn_subsample(tempBehav,tempBefIdx,nRep,nProbe);
                    probeByDayBin.aft(i,1) = nanmean(tempAftAcc); probeByDayBin.aft(i,2) = nanmean(abs(tempAftBias));
                end 

            end
            function [tempBias,tempAcc] = fn_subsampleBlock(behav,selIdx,nSamp)
                tempBias = nan(1,length(selIdx)); tempAcc= nan(1,length(selIdx));
                if nSamp>0
                    tempStartSelFlag = zeros(nSamp,length(selIdx));
                    for m = 1:nSamp
                        selIdxShifted = selIdx+m-1; tempStartSelFlag(m,:) =sum(selIdx==selIdxShifted',1);

                    end
                    tempStartSelFlag = sum(tempStartSelFlag,1); startSelIdx = find(tempStartSelFlag==nSamp);
                    
                    for n = 1:length(startSelIdx)
                        tempIdx = startSelIdx:startSelIdx+nSamp-1;
                        [tempBias(n), tempAcc(n), ~, ~, ~] = fn_getAccBias(behav.stimulus(tempIdx), ...
                            behav.responseType(tempIdx)==1,behav.responseType(tempIdx)==0);
                    end
                end
            end
            
        end
        % REINF VS. PROBE LICK COMPARISON, (BEF VS. AFT) OR (SUBSAMPLE WITHIN DAY)
        function [nLickReinf, nLickProbe, lickTimeReinf, lickTimeProbe,tempLickTimeAllReinf,tempLickTimeAllProbe] = computeLickAlignByTrial(obj,compType,accThre,windowSize)
            if nargin == 1; accThre = [0.7 nan]; windowSize= 400; end
            [~, learningOnset] = obj.fn_findLearningOnset(accThre,windowSize);
            selTrial = learningOnset-300:learningOnset+300;
            nTrials = size(obj.behav,1); 
            selTrial(selTrial<=0) = []; selTrial(selTrial>nTrials) = [];
            dayStart = obj.behav.day(selTrial(1)); dayEnd = obj.behav.day(selTrial(end));

            if strcmp(compType,'befAft')
                [nLickReinf, nLickProbe, lickTimeReinf, lickTimeProbe,tempLickTimeAllReinf,tempLickTimeAllProbe] = fn_lickByDay(obj,[dayStart dayEnd],compType);
            elseif strcmp(compType,'day')
                [nLickReinf, nLickProbe, lickTimeReinf, lickTimeProbe,tempLickTimeAllReinf,tempLickTimeAllProbe] = fn_lickByDay(obj,[dayStart dayEnd],compType);
            end
        end
        function [lickReinf, lickProbe, lickTimeReinf, lickTimeProbe,tempLickTimeAllReinf,tempLickTimeAllProbe] = fn_lickByDay(obj,selDay,compType)
            nLickReinf = []; nLickProbe = [];
            dayStart = selDay(1); dayEnd = selDay(2);nDay = length(dayStart:dayEnd);
            trialLim(1) = find(obj.behav.day==dayStart,1); trialLim(2) = find(obj.behav.day==dayEnd,1,'last'); %nBefAft = 30;

            probeByDay.probe = nan(nDay,2); probeByDay.reinf = nan(nDay,2);
            probeByDay.bef = nan(nDay,2); probeByDay.aft = nan(nDay,2); 
            for i = 1:nDay
                tempDay = dayStart+i-1; 
                tempBehav = obj.behav(obj.behav.day == tempDay,:); 
                probeFlag = tempBehav.goodProbe;
                corrFlag = tempBehav.responseType==1;
                nLickReinf(i,1) = nanmean(tempBehav.nlick(~probeFlag & ~tempBehav.badReinf & corrFlag));
                nLickReinf(i,2) = nanmean(tempBehav.nlick(~probeFlag & ~tempBehav.badReinf & ~corrFlag));
                nLickProbe(i,1) = nanmean(tempBehav.nlick(probeFlag & corrFlag));
                nLickProbe(i,2) = nanmean(tempBehav.nlick(probeFlag & ~corrFlag));

                reinfCorr = ~probeFlag & ~tempBehav.badReinf & corrFlag;
                reinfInCorr = ~probeFlag & ~tempBehav.badReinf & ~corrFlag;

                [lickReinf{1}, lickTimeReinf{1},tempLickTimeAllReinf{1}] = getLickTime(tempBehav,reinfCorr);
                [lickReinf{2}, lickTimeReinf{2},tempLickTimeAllReinf{2}] = getLickTime(tempBehav,reinfInCorr);

                probeCorr = probeFlag & corrFlag;
                probeInCorr = probeFlag & ~corrFlag;

                [lickProbe{1}, lickTimeProbe{1},tempLickTimeAllProbe{1}] = getLickTime(tempBehav,probeCorr);
                [lickProbe{2}, lickTimeProbe{2},tempLickTimeAllProbe{2}] = getLickTime(tempBehav,probeInCorr);
                
                
            end

            function [tempLickCount,tempLickTime,tempLickTimeAll] = getLickTime(behav,trialFlag)
                tempStr = behav.lickTimeStr(trialFlag); tempLickCount = []; tempLickTime = [];tempLickTimeAll = {};
                tempRT = behav.reactionTime(trialFlag);
                for j = 1:length(tempStr)
                    temp = strsplit(tempStr{j}); 
                    tempTime = cellfun(@(x)(str2double(x)),temp,'UniformOutput',false);
                    if iscell(tempTime); tempTime = cell2mat(tempTime); end
                    selLickIdx = (tempTime-tempRT(j))<1 & (tempTime-tempRT(j))>0;
                    tempLickCount(j) = sum(selLickIdx); 
                    temp = min(tempTime(selLickIdx)); 
                    if isempty(temp); temp = nan; end 
                    tempLickTime(j) = temp; 
                    tempLickTimeAll{j} = tempTime-tempRT(j); 
                end


            end
        end

        % FUNCTION -- find the onset of learning
        function [trialBin, learningOnset] = fn_findLearningOnset(obj,accThre, binSize)
            reinfAccThre = accThre(1); probeAccThre = accThre(2);  probeBin = obj.binProbe(true);
            tempProbeBlockStart = find(probeBin.probeAcc > probeAccThre,1);
            if isempty(tempProbeBlockStart) % if no probe acc threshold is given
                probeStartIdx = nan; 
            elseif tempProbeBlockStart > 1 % if later probe bins reach threhold, take the mid point with previous bin
                probeStartIdx = round((probeBin.onIdx(tempProbeBlockStart) + probeBin.offIdx(tempProbeBlockStart-1))/2);
            else  % if first probe bin reach threshold, then use the beginning of this bin
                probeStartIdx = probeBin.onIdx(tempProbeBlockStart);
            end 
            reinfAccBin = 300;
            [~,reinfAcc] = fn_getAccBiasSmooth(obj.behav.stimulus,obj.behav.responseType,reinfAccBin);
            reinfStartIdx = find(reinfAcc > reinfAccThre,1);
            if isempty(reinfStartIdx) % if no reinf acc threshold is given
                reinfStartIdx = nan; end 
            % determine the learning onset and give the bins
            learningOnset = min(reinfStartIdx,probeStartIdx);
            %learningOnset = find(reinfAcc > 0.6,1);
            nTrials = size(obj.behav,1);
            trialBin1 = fliplr(learningOnset:-binSize:1); 
            if length(trialBin1)> 1; trialBin1(1) = 1; else; trialBin1 = [1 learningOnset]; end
            trialBin2 = learningOnset:binSize:nTrials; trialBin2(end) = nTrials;
            trialBin = [trialBin1 trialBin2(2:end)];
        end       
        % FUNCTION -- get data and correct for different things
        function obj = removeProbe(obj)
            missFlag = (obj.behav.goodProbe == 1); obj.behav = obj.behav(~missFlag,:);
        end      
        function obj = correctWheelTrial(obj)
            if strcmp(obj.ops.inputType,'cell')
                for i = 1:length(obj.stimulus)
                    if ~isempty(obj.wheel.wheelSoundOn{i})
                        obj.wheel.wheelSoundOn{i} = nansum(~isnan(obj.wheel.wheelSoundOn{i}),2);
                    else
                        obj.wheel.wheelSoundOn{i}=nan(length(obj.stimulus{i}),1);
                    end

                    if ~obj.wheel.wheelSoundOnCheckFlag{i}
                        obj.wheel.wheelPreSound{i} = nan(length(obj.stimulus{i}),1);
                        obj.wheel.wheelSoundOn{i} = nan(length(obj.stimulus{i}),1);
                    end
                end
            else; msgbox(['Multiday type is ' obj.ops.multidayType ', need to be cell'], 'ERROR MESSAGE');
            end
        end
        function obj = correctMultiTaskStim(obj) 
            if isfield(obj.ops,'correctionMultiTaskName') && ~isempty(obj.ops.correctionMultiTaskName)
                correctionTaskFlag = contains(obj.behav.trainingType,obj.ops.correctionMultiTaskName);
                temp = obj.behav.stimulus(correctionTaskFlag);
                temp1Flag = (temp ==1 | temp ==2) ; temp2Flag = (temp ==3 | temp ==4) ; 
                temp(temp1Flag) = temp(temp1Flag) + 2; temp(temp2Flag) = temp(temp2Flag) - 2;
                obj.behav.stimulus(correctionTaskFlag) = temp;
            end
        end
        function obj = getPsyTrack(obj)

            if isfield(obj.ops,'psytrackPath')
                model = load([obj.ops.psytrackPath filesep obj.ops.mouse 'psytrack_SBARpRn_nPrev1.mat']);
                prevBehav = load(['C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData\trialData' filesep obj.ops.mouse '_nPrev5.mat']);
                modelAhW = model.wMode(1,:)'; modelRnW = model.wMode(2,:)'; modelRpW = model.wMode(3,:)';
                modelBiasW = model.wMode(4,:)';
                modelAccW = model.wMode(5,:)';
                
                trialDiff = size(obj.behav,1) - size(model.wMode,2);
                model = attachTrials(model,trialDiff,{'wMode'},false);
                prevBehav = attachTrials(prevBehav,trialDiff,{'stimulus','actionH','actionXnegRewardH','actionXposRewardH'},true);
                %{
                if trialDiff>0
                    modelBiasW = [repmat(modelBiasW(1),trialDiff,1); modelBiasW];
                    modelAccW = [repmat(modelAccW(1),trialDiff,1); modelAccW];
                else
                    modelBiasW = modelBiasW(1:size(obj.behav,1));
                    modelAccW = modelAccW(1:size(obj.behav,1));
                end
                %}
                obj.behav.modelBias = 1-fn_logistic(model.wMode(:,4))*2; % bias given by the model by weight
                obj.behav.modelAcc = fn_logistic(model.wMode(:,5)); % acccuracy on the task, given by stimulus weight only
                
                obj.behav.modelPred = fn_logistic( prevBehav.actionH(:,1) .* model.wMode(:,1) + ...
                    prevBehav.actionXnegRewardH(:,1) .* model.wMode(:,2)...
                    + prevBehav.actionXposRewardH(:,1) .* model.wMode(:,3) + model.wMode(:,4)...
                    + prevBehav.stimulus .* model.wMode(:,5) ); % predicted choice with probability
                obj.behav.modelPredCorr = (1-obj.behav.modelPred) .* (prevBehav.stimulus == -1) +...
                    obj.behav.modelPred .* (prevBehav.stimulus == 1); % probability of model choosing the correct choice
                [bias,acc_L,acc_R] = fn_getBiasModel(prevBehav.stimulus , obj.behav.modelPredCorr,obj.ops.learningCurveBin);
                obj.behav.modelPredAcc = (acc_L+acc_R)/2; % probability of model choosing the correct choice, smoothed and corrected for stim prob
                obj.behav.modelPredBias= bias;

                obj.behav.modelPred_removeA = fn_logistic( prevBehav.actionH(:,1) .* model.wMode(:,1) + ...
                    prevBehav.actionXnegRewardH(:,1) .* model.wMode(:,2)...
                    + prevBehav.actionXposRewardH(:,1) .* model.wMode(:,3) ...
                    + prevBehav.stimulus .* model.wMode(:,5) );
                obj.behav.modelPred_removeACorr = (1-obj.behav.modelPred_removeA) .* (prevBehav.stimulus == -1) +...
                    obj.behav.modelPred_removeA .* (prevBehav.stimulus == 1);
                [~,acc_L,acc_R] = fn_getBiasModel(prevBehav.stimulus ,obj.behav.modelPred_removeACorr,obj.ops.learningCurveBin);
                obj.behav.modelPred_removeAAcc = (acc_L+acc_R)/2;
                
                obj.behav.modelPred_removeH = fn_logistic( model.wMode(:,4)...
                    + prevBehav.stimulus .* model.wMode(:,5) );
                obj.behav.modelPred_removeHCorr = (1-obj.behav.modelPred_removeH) .* (prevBehav.stimulus == -1) +...
                    obj.behav.modelPred_removeH .* (prevBehav.stimulus == 1);
                [~,acc_L,acc_R] = fn_getBiasModel(prevBehav.stimulus ,obj.behav.modelPred_removeHCorr,obj.ops.learningCurveBin);
                obj.behav.modelPred_removeHAcc = (acc_L+acc_R)/2;
                
                
                obj.behav.modelPred_removeAll = fn_logistic( prevBehav.stimulus .* model.wMode(:,5) );
                obj.behav.modelPred_removeAllCorr = (1-obj.behav.modelPred_removeAll) .* (prevBehav.stimulus == -1) +...
                    obj.behav.modelPred_removeAll .* (prevBehav.stimulus == 1);
                [~,acc_L,acc_R] = fn_getBiasModel(prevBehav.stimulus ,obj.behav.modelPred_removeAllCorr,obj.ops.learningCurveBin);
                obj.behav.modelPred_removeAllAcc = (acc_L+acc_R)/2;
                
            end
            
            function S = attachTrials(S,trialDiff,targetFields,nanflag)
            
                for i = 1:length(targetFields)
                    if size(S.(targetFields{i}),1) < size(S.(targetFields{i}),2)
                        S.(targetFields{i}) = S.(targetFields{i})';
                    end
                    if trialDiff>0
                        if ~nanflag
                            S.(targetFields{i}) = cat(1,repmat(S.(targetFields{i})(1,:),trialDiff,1), S.(targetFields{i}));
                        else
                            S.(targetFields{i}) = cat(1,repmat([nan],trialDiff,size(S.(targetFields{i}),2)), S.(targetFields{i}));
                        end
                    else
                        S.(targetFields{i}) = S.(targetFields{i})(1:size(obj.behav,1),:);
                    end

                end
            end
        end
        % BIASED BLOCK COMPUTATION
        function obj = getBiasBlock_oldNoCorrection(obj)
            actL = obj.action == 1;biasBin = 20;
            actL = smoothdata(actL,'movmean',biasBin);
            actAxis = biasBin/2+1:length(actL)-biasBin/2;

            biasThreshold = 0.15; trialBlockThreshold = 5;
            biasL = find(actL >= (biasThreshold+0.5) ); biasR = find(actL <= (0.5-biasThreshold));
            obj.biasBlockL = selectBlock(biasL); obj.biasBlockR = selectBlock(biasR);
            
            function biasBlock = selectBlock(bias)
                bias_incre = diff(bias); temp = bias_incre(2:end);
                startFlag = find(bias_incre>=5); 
                biasBlock.start = bias([1; startFlag+1]) ; 
                biasBlock.end = bias([startFlag; length(bias)]) ;
                biasBlock.len= biasBlock.end - biasBlock.start+1;

                threFlag = biasBlock.len>=trialBlockThreshold;
                biasBlock = structfun(@(x)(x(threFlag)),biasBlock,'UniformOutput', false);
            end
            
        end
        function obj = getBiasBlock(obj,actionBiasThre)
            %actionBiasThre = 0.4; stimCorrThre = 0.50; transTrialThreshold = 50;
            trialBin = 30;
            %bias = obj.behav.modelBias; 
            bias = obj.behav.bias;
            
                
            stimulus = obj.behav.stimulus;
            stimL = (stimulus == 1); stimL = smoothdata(stimL,'movmean',trialBin);
            [obj.biasBlock.stateFlag, obj.biasBlock.blockL, obj.biasBlock.blockR, obj.biasBlock.blockU]=...
                fn_detectBlock(bias,0, actionBiasThre,'blockLenThre',10,'blockIntervalThre',10);
            
            [obj.biasBlock.trans, obj.biasBlock.transID] = fn_getTransition(obj.biasBlock);
            %thresholdFlag = (L2U(:,3)-L2U(:,2)) <= transTrialThreshold; L2U = L2U(thresholdFlag,:);
            %thresholdFlag = (U2L(:,3)-U2L(:,2)) <= transTrialThreshold; U2L = U2L(thresholdFlag,:);
            %obj.biasBlock.L2R = stimProbCorrection(L2R,stimL, stimCorrThre); 
            %obj.biasBlock.R2L = stimProbCorrection(R2L,stimL, stimCorrThre);

            function L2R = stimProbCorrection(L2R,stimL, stimThre)
                if ~isempty(L2R)
                    for i = 1:size(L2R,1)
                        tempStim = stimL(L2R(i,1):L2R(i,4));
                        tempMin = min(tempStim); tempMax = max(tempStim);
                        tempFlag(i) = (tempMin<(0.5-stimThre)) | (tempMax > (0.5+stimThre));  
                    end
                    L2R = L2R (~tempFlag,:);
                end
            end
        end     
        function obj = getBiasBlock2(obj)

            mouseRand = rand(1,1000000); mouseRand(mouseRand<0.5) = -1; mouseRand(mouseRand>0.5) = 1;
            biasRand = smoothdata(mouseRand,'movmean',obj.ops.learningCurveBin);
            [~, ~, ~,bareaL] = fn_getBlockOnOff(biasRand>0,biasRand);
            [~, ~, ~,bareaR] = fn_getBlockOnOff(biasRand<0,biasRand);
            barea = [bareaL bareaR]; areaThre = prctile(barea,95);
            [onL, offL, idxL,areaL] = fn_getBlockOnOff(obj.behav.bias>0,obj.behav.bias);
            [onR, offR, idxR,areaR] = fn_getBlockOnOff(obj.behav.bias<0,obj.behav.bias);
            selFlagL = find(areaL > prctile(barea,95));
            selFlagR = find(areaR > prctile(barea,95));
            
            obj.biasBlock.blockL.start = onL(selFlagL); obj.biasBlock.blockL.end = offL(selFlagL);
            idxL = idxL(selFlagL); obj.biasBlock.blockL.len = cellfun(@length,idxL,'UniformOutput',true);
            obj.biasBlock.blockR.start = onR(selFlagR); obj.biasBlock.blockR.end = offR(selFlagR);
            idxR = idxR(selFlagR); obj.biasBlock.blockR.len = cellfun(@length,idxR,'UniformOutput',true);
            
            % calculate the stateFlag for the 1st time
            obj.biasBlock.stateFlag = zeros(1,size(obj.behav,1)); 
            for i = 1:length(idxL); obj.biasBlock.stateFlag(idxL{i}) = 1; end
            for i = 1:length(idxR); obj.biasBlock.stateFlag(idxR{i}) = -1; end

            [obj.biasBlock.blockU.start, obj.biasBlock.blockU.end, idxU,~] = fn_getBlockOnOff(obj.biasBlock.stateFlag==0,obj.biasBlock.stateFlag);
            obj.biasBlock.blockU.len =  cellfun(@length,idxU,'UniformOutput',true);
            % FIX SHORT UNBIASED BLOCKS AND BIASED BLOCKS WITH ALL NANS
            uBlockThre = 10; tempB = [];
            
            for i = 1:length(obj.biasBlock.blockU.start)
                tempB(i) = sum(isnan(obj.behav.bias(obj.biasBlock.blockU.start(i):obj.biasBlock.blockU.end(i)))); 
                tempBTotal(i) = length(obj.behav.bias(obj.biasBlock.blockU.start(i):obj.biasBlock.blockU.end(i))); 
            end
            tempNanFlag = tempB>1; 
            uBlockFlag = find(obj.biasBlock.blockU.len < uBlockThre | (tempNanFlag & tempB==tempBTotal));
            
            for i = 1:length(uBlockFlag)
                tempOn = obj.biasBlock.blockU.start(uBlockFlag(i))-1; if tempOn<1; tempOn = 1; end
                tempOff = obj.biasBlock.blockU.end(uBlockFlag(i))+1; if tempOff>length(obj.behav.bias); tempOff=length(obj.behav.bias); end
                tempBiasDir = sum(obj.behav.bias(obj.biasBlock.blockU.start(uBlockFlag(i)):obj.biasBlock.blockU.end(uBlockFlag(i))));
                
                tempBefL = find(obj.biasBlock.blockL.end == tempOn);
                tempAftL = find(obj.biasBlock.blockL.start == tempOff);
                tempBefR = find(obj.biasBlock.blockR.end == tempOn);
                tempAftR = find(obj.biasBlock.blockR.start == tempOff);
                if  ~isempty(tempBefL); concatBef = 'L'; elseif ~isempty(tempBefR); concatBef = 'R'; else; concatBef = ''; end
                if  ~isempty(tempAftL); concatAft = 'L'; elseif ~isempty(tempAftR); concatAft = 'R'; else; concatBef = ''; end
                if isnan(tempBiasDir) || tempBiasDir==0; biasFlag = concatAft; elseif tempBiasDir > 0;  biasFlag = 'L'; else biasFlag = 'R'; end

                if strcmp(biasFlag,'L') && (strcmp(concatBef,'L') || strcmp(concatAft,'L')); concatFlag = 'L';
                elseif strcmp(biasFlag,'R') && (strcmp(concatBef,'L') && strcmp(concatAft,'L'));concatFlag = 'L';
                else concatFlag = 'R';
                end

                if strcmp(concatFlag,'L') % left biased block
                    tempBef = find(obj.biasBlock.blockL.end == tempOn);
                    tempAft = find(obj.biasBlock.blockL.start == tempOff);
                    if ~isempty(tempBef)
                        obj.biasBlock.blockL.end(tempBef) = tempOff;
                        obj.biasBlock.blockL.len(tempBef) = tempOff-obj.biasBlock.blockL.start(tempBef)+1;
                    elseif ~isempty(tempAft)
                        obj.biasBlock.blockL.start(tempAft) = tempOn;
                        obj.biasBlock.blockL.len(tempAft) = obj.biasBlock.blockL.end(tempAft) - tempOn +1;
                    end
                elseif strcmp(concatFlag,'R') % right biased block
                    tempBef = find(obj.biasBlock.blockR.end == tempOn);
                    tempAft = find(obj.biasBlock.blockR.start == tempOff);
                    if ~isempty(tempBef)
                        obj.biasBlock.blockR.end(tempBef) = tempOff;
                        obj.biasBlock.blockR.len(tempBef) = tempOff-obj.biasBlock.blockR.start(tempBef)+1;
                    elseif ~isempty(tempAft)
                        obj.biasBlock.blockR.start(tempAft) = tempOn;
                        obj.biasBlock.blockR.len(tempAft) = obj.biasBlock.blockR.end(tempAft) - tempOn +1;

                    end
                end
            end
            
            i = 2; 
            while i <= length(obj.biasBlock.blockL.start)
                if any(obj.biasBlock.blockL.start(i) == obj.biasBlock.blockL.end(i-1))
                    obj.biasBlock.blockL.start(i) = [];
                    obj.biasBlock.blockL.end(i-1) = [];
                    obj.biasBlock.blockL.len(i) = [];
                    obj.biasBlock.blockL.len(i-1) = obj.biasBlock.blockL.end(i-1)-obj.biasBlock.blockL.start(i-1)+1;
                else 
                    i = i+1;
                end
                
            end

            i = 2; 
            while i <= length(obj.biasBlock.blockR.start)
                if any(obj.biasBlock.blockR.start(i) == obj.biasBlock.blockR.end(i-1))
                    obj.biasBlock.blockR.start(i) = [];
                    obj.biasBlock.blockR.end(i-1) = [];
                    obj.biasBlock.blockR.len(i) = [];
                    obj.biasBlock.blockR.len(i-1) = obj.biasBlock.blockR.end(i-1)-obj.biasBlock.blockR.start(i-1)+1;
                else
                    i = i+1;
                end
                
            end
            % Re-calculate the unbiased blocks
            obj.biasBlock.stateFlag = zeros(1,size(obj.behav,1)); 
            for i = 1:length(obj.biasBlock.blockL.start); obj.biasBlock.stateFlag(obj.biasBlock.blockL.start(i):obj.biasBlock.blockL.end(i)) = 1; end
            for i = 1:length(obj.biasBlock.blockR.start); obj.biasBlock.stateFlag(obj.biasBlock.blockR.start(i):obj.biasBlock.blockR.end(i)) = -1; end
            [obj.biasBlock.blockU.start, obj.biasBlock.blockU.end, idxU,~] = fn_getBlockOnOff(obj.biasBlock.stateFlag==0,obj.biasBlock.stateFlag);
            obj.biasBlock.blockU.len =  cellfun(@length,idxU,'UniformOutput',true);
            % GET BIAS TRANSITION IDENTITY
            [obj.biasBlock.trans, obj.biasBlock.transID] = fn_getTransition(obj.biasBlock);

            % CALCULATE THE MEAN BIAS WITHIN EACH BLOCK
            [obj.biasBlock.blockL.acc,obj.biasBlock.blockL.bias,obj.biasBlock.blockL.trial,obj.biasBlock.blockL.peakIdx] = getBlockData(obj.biasBlock.blockL,obj.behav);
            [obj.biasBlock.blockR.acc,obj.biasBlock.blockR.bias,obj.biasBlock.blockR.trial,obj.biasBlock.blockR.peakIdx] = getBlockData(obj.biasBlock.blockR,obj.behav);
            [obj.biasBlock.blockU.acc,obj.biasBlock.blockU.bias,obj.biasBlock.blockU.trial,obj.biasBlock.blockU.peakIdx] = getBlockData(obj.biasBlock.blockU,obj.behav);

            disp('Done')
            function [acc,bias,trial,peakIdx] = getBlockData(block,behav)
                trial = block.start/2 + block.end/2;
                acc = []; bias = []; peakIdx = {};
                for k = 1:length(block.start)
                    tempOnOff = floor(block.len(k)/3); tempOnOff = min([tempOnOff 20]);
                    tempIdx = block.start(k)+tempOnOff : block.end(k)-tempOnOff;
                    acc(k) = nanmean(behav.acc(tempIdx));
                    bias(k) = nanmean(abs(behav.bias(tempIdx)));
                    peakIdx{k} = tempIdx;
                end
            end

        end        
    end
end

function [bias, acc] = fn_calculateBiasBySample(mat,nSample,nProbe)
    tempSample = zeros(nSample,6); 
    tempSample(:,1:2) = fn_generateSampleByCount(mat(1:2),nSample,nProbe(1));
    tempSample(:,4:5) = fn_generateSampleByCount(mat(4:5),nSample,nProbe(2));
    [bias,acc,~] = fn_getAccBiasByCount(tempSample); 
end

function [sample] = fn_generateSampleByCount(mat,nSample,nProbe)
    tempCount = sum(mat); tempSample = 1:tempCount;
    tempSample(1:mat(1)) = 1; tempSample(mat(1)+1:end) = -1;
    sample = zeros(nSample,2);
    for m = 1:nSample
        try
            tempRand = randperm(length(tempSample),nProbe);
        catch 
            tempRand = 1:length(tempSample); % if the number of trials is smaller than probe, do not permute
        end
        sample(m,1) = sum(tempSample(tempRand) == 1);
        sample(m,2) = sum(tempSample(tempRand) == -1);
    end

end

function [tempBias,tempAcc,tempDP, tempCri] = fn_subsample(behav,selIdx,nRep,nSamp)
    tempBias = nan(1,nRep); tempAcc= nan(1,nRep); tempDP= nan(1,nRep); tempCri= nan(1,nRep);
    if nSamp>0
        for n = 1:nRep
            tempRand1 = randperm(length(selIdx),nSamp); tempRandIdx = selIdx(tempRand1);
            [tempBias(n), tempAcc(n), tempDP(n), tempCri(n), ~] = fn_getAccBias(behav.stimulus(tempRandIdx), ...
                behav.responseType(tempRandIdx)==1,behav.responseType(tempRandIdx)==0);
        end
    end
end
