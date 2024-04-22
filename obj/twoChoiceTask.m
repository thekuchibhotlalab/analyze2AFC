classdef twoChoiceTask
    %TASKOBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % CONVENTIONS -- 
        % Stim in 1,2,3,4... (or float nnmbers)
        % Context in 1,2,3...
        % ResponsType in 0(miss), 1(correct), 2(incorrect)
        % Action - 0(no action), 1(choice 1), 2(choice 2)
        behav
        behavVar
        ops
        
    end
    
    methods
        %-----------------------------CONSTRUCTOR METHODS-----------------------------------------
        function [obj,inStruct] = twoChoiceTask(inStruct,varargin)
            %TASKOBJ Construct an instance of this class
            %   Detailed explanation goes here       
            
            p = fn_inputParser(); p.parse(varargin{:});
            inputNames = fieldnames(inStruct);
            taskObjPropName = properties('twoChoiceTask');
            % ---------- 1.1 Write variables into object ---------
            %addpath([p.Results.opsPath filesep]); feval([p.Results.mouse '_ops']); rmpath([p.Results.opsPath filesep]);
            
            obj.ops = p.Results.ops;
            obj.ops.mouse = p.Results.mouse;
            obj.ops.verbose = p.Results.verbose;
            if ~isfield(obj.ops,'accInterpolation'); obj.ops.accInterpolation = false; end
            if ~isfield(obj.ops,'learningCurveBin') || isempty(obj.ops.learningCurveBin); obj.ops.learningCurveBin = 100; end
            
            if ~isfield(obj.ops,'preprocessedInput') || ~obj.ops.preprocessedInput; inStruct = fn_preprocessInput(inStruct,obj.ops); end
            
            obj.behav = fn_struct2table(inStruct,'fieldKeys',getBehTableKeys(obj.ops.behavType));
            obj.behavVar = struct();
            obj.behav.actionRate = smoothdata(obj.behav.action~=0,'movmean',obj.ops.learningCurveBin);
            obj.behav.masterTrialnum = (1:length(obj.behav.stimulus))';
            tempTrialnumSession = obj.behav.masterTrialnum;
            for i = 1:obj.behav.day(end)
                tempDayFlag = find(obj.behav.day==i);
                tempTrialnumSession(tempDayFlag) = tempTrialnumSession(tempDayFlag) - tempDayFlag(1)+1;
            end
            obj.behav.masterTrialnumSession = tempTrialnumSession;


            obj = getProbe(obj);

            function p = fn_inputParser()
                p = inputParser;
                p.KeepUnmatched = true;
                arg = {'ops','mouse','verbose';...
                        []     , []    ,false};
                cellfun(@(x,y)(p.addParameter(x,y)),arg(1,:),arg(2,:));
            end    

            function obj = getProbe(obj)
                obj.behav.probe = (obj.behav.context == 3);
                obj = correctForBadProbe(obj);
            end
            
        end
        %-----------------------------METHODS SPECIFIC TO TWOCHOICE TASK-----------------------------------------
        
        function obj = correctForBadProbe(obj)
            T = separateDay(obj);
            for i = 1:length(T)
                nProbe = sum(T{i}.probe);
                if mod(nProbe,10) == 0; T{i}.goodProbe = T{i}.probe;
                else; T{i}.goodProbe = zeros(size(T{i}.probe));
                end 
                [probeOn, probeOff] = fn_getBlockOnOff(T{i}.goodProbe,T{i}.goodProbe);
                probeIdx = find(T{i}.goodProbe); probeOnIdx = find(probeOn); probeOffIdx = find(probeOff);
                if isempty(probeOnIdx); reinfBefIdx = [];reinfAftIdx = [];
                elseif (probeOffIdx(1) - probeOnIdx(1)) == 9
                    reinfBefIdx = probeIdx - 10; reinfAftIdx = probeIdx + 10;
                elseif (probeOffIdx(1) - probeOnIdx(1)) == 19
                    reinfBefIdx = probeIdx - 30; reinfAftIdx = probeIdx + 30;
                else; reinfBefIdx = [];reinfAftIdx = [];
                end
                reinfBefIdx(reinfBefIdx < 1) = [];
                reinfAftIdx(reinfAftIdx > length(T{i}.probe)) = [];
                T{i}.reinfBef = zeros(size(T{i}.probe)); T{i}.reinfBef(reinfBefIdx) = 1;
                T{i}.reinfAft = zeros(size(T{i}.probe)); T{i}.reinfAft(reinfAftIdx) = 1; 
            end
            temp = vertcat(T{:});
            obj.behav = temp;
        end
        
        function T = separateDay(obj)
            days = unique(obj.behav.day);
            T = cellfun(@(x)(obj.behav(obj.behav.day==x,:)),num2cell(days),'UniformOutput',false);
        end
        
        function obj = removeMiss(obj)
            missFlag = (obj.behav.action == 0); obj.behav = obj.behav(~missFlag,:);
            obj.behavVar = fn_structfun(@(x)(x(~missFlag,:)),obj.behavVar);
            obj.behav.trialnum = (1:length(obj.behav.stimulus))';
            tempTrialnumSession = obj.behav.trialnum;
            for i = 1:obj.behav.day(end)
                tempDayFlag = find(obj.behav.day==i);
                tempTrialnumSession(tempDayFlag) = tempTrialnumSession(tempDayFlag) - tempDayFlag(1)+1;
            end
            obj.behav.trialnumSession = tempTrialnumSession;
        end
        
        %function obj = selectDay(obj,dayFlag)
        %    if strcmp(obj.ops.inputType,'cell')
        %        obj = fn_readStructBySelection(obj,'selectFlag',dayFlag);
        %        obj.dayLen.total = obj.dayLen.total(dayFlag);
        %        obj.dayLen.noMiss = obj.dayLen.noMiss(dayFlag);
        %    else; msgbox(['Input type is ' obj.ops.inputType ', need to be cell'], 'ERROR MESSAGE');
        %    end
        %end
        
        %function obj = concatenateDay(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
        %    if strcmp(obj.ops.inputType,'cell')
        %        propCell = properties('twoChoiceTask');
        %        for i = 1:length(propCell)
        %            if iscell(obj.(propCell{i})); obj.(propCell{i}) = fn_cell2mat(obj.(propCell{i}),1);end
        %        end
        %    else; msgbox(['Input type is ' obj.ops.inputType ', need to be cell'], 'ERROR MESSAGE');
        %    end
        %end

%         function obj = removeIdx(obj,idx)
%             if strcmp(obj.ops.multidayType,'mat')       
%                 obj = fn_objfun(@(x)fn_removeIdx(x,idx),obj,'verbose',obj.ops.verbose,'funcName', 'REMOVE_IDX');
%             else
%                 msgbox(['Multiday type is ' obj.ops.multidayType ', need to be mat'], 'ERROR MESSAGE');
%             end 
%         end


        function obj = getAcc(obj,learningCurveBin)
            if nargin == 1; learningCurveBin = obj.ops.learningCurveBin; end

            [obj.behav.bias,obj.behav.acc,obj.behav.dprime,obj.behav.crit,obj.behav.acc1,obj.behav.acc2] = ...
                fn_getAccBiasSmooth(obj.behav.stimulus,obj.behav.responseType,learningCurveBin);

            tempNanFlag1 = find(isnan(obj.behav.acc1) & ~isnan(obj.behav.acc2));
            if ~isempty(tempNanFlag1)
                if obj.ops.accInterpolation
                    obj.behav.acc1(tempNanFlag1) = interp1(1:length(obj.behav.acc1),obj.behav.acc1,(tempNanFlag1));
                    obj.behav.acc1(obj.behav.acc1>1) = 1; obj.behav.acc1(obj.behav.acc1<0) =0;
                    disp('acc1 interpolation happen!')
                else; disp('acc1 nan detected!');end
            end
            tempNanFlag2 = find(~isnan(obj.behav.acc1) & isnan(obj.behav.acc2));
            if ~isempty(tempNanFlag2) 
                if obj.ops.accInterpolation
                    obj.behav.acc2(tempNanFlag2) = interp1(1:length(obj.behav.acc2),obj.behav.acc2,(tempNanFlag2));
                    obj.behav.acc2(obj.behav.acc2>1) = 1; obj.behav.acc2(obj.behav.acc2<0) =0;
                    disp('acc2 interpolation happen!')
                else; disp('acc2 nan detected!');end
            end
            if ~isempty(tempNanFlag1) || ~isempty(tempNanFlag2) % recalculate the 
                obj.behav.bias = obj.behav.acc1 - obj.behav.acc2;
                obj.behav.acc = obj.behav.acc1/2 + obj.behav.acc2/2;
                
            end
        end

        function [taskAcc,taskBias,taskDP,taskAcc1,taskAcc2] = getAccMultiTask(obj)

            nTask = max(unique(obj.behav.stimulus))/2;taskAcc = {}; taskAcc1 = {}; taskAcc2 = {}; taskBias = {};taskDP = {};
            for i = 1:nTask
                tempStimL = 1 + (i-1)*2; tempStimR = 2 + (i-1)*2;
                taskFlag = (obj.behav.stimulus == tempStimL | obj.behav.stimulus == tempStimR);

                tempStim = obj.behav.stimulus; tempStim(~taskFlag) = nan;
                tempResp = obj.behav.responseType; tempResp(~taskFlag) = nan;
                [taskBias{i},taskAcc{i},taskDP{i},~,taskAcc1{i},taskAcc2{i}] = fn_getAccBiasSmooth(tempStim,tempResp,obj.ops.learningCurveBin);
                %taskAcc{i} = (acc_L +acc_R)/2;
                
                
                %dprimeAccThre = 0.95;
                %acc_L(acc_L>dprimeAccThre) = dprimeAccThre; acc_R(acc_R>dprimeAccThre) = dprimeAccThre;
                %taskDP{i} = norminv(acc_L,0,1) - norminv(1-acc_R,0,1);
            end
                        
            
        end
        

    end
end

