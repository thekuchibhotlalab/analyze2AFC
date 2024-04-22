function  readCoulbournData(mousename)
global sep mouse;
sep = '\';
mouse = mousename;
dataPath = 'C:\Users\zzhu34\Documents\tempdata\octoData';
writeDataPath = [dataPath sep 'trialData' sep]; mkdir(writeDataPath);
wheelTxtPath = ['C:\Users\zzhu34\Documents\tempdata\octoData\CoolTermData' sep mouse sep];

readDataPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\coulbournData';
mousePath = [readDataPath sep mouse sep];
%filenames = dir([mousePath '*.csv']);

diary([ writeDataPath sep 'log.txt'])
%% read animal name and date for each file
dirName = fn_readDir(mousePath);
trialData = {};
for i = 1:length(dirName)
   trialData{i} = readFiles([mousePath sep dirName{i}],dirName{i},wheelTxtPath); 
end
if length(trialData) > 1; trialData = fn_catStructField(2, trialData{:}); 
else; trialData = trialData{1};end
trialData = fn_sortStructByFieldKey(trialData,'date');

save([writeDataPath sep mouse '.mat'],'trialData');
diary off;
end
%% All the functions

function trialData = readFiles(csvPath,dirName,wheelTxtPath)
global sep mouse; 
filenames = dir([csvPath sep '*.csv']);
stateKey = fn_readCoulbournStateKey();
trialName = []; responseName = [];

trialData = struct();
trainDay = 0; tempDate = [];
tempEntryIdx = []; tempEntryName = []; tempEntryTime = []; 
tempInputIdx = []; tempInputName = []; tempInputTime = []; 
tempVar4Idx = []; tempVar4Name = []; tempVar4Time = []; 
expDate = {}; expDateSimp = {};
for i = 1:length(filenames)
    filename = filenames(i).name;
    filenameSplit = strsplit(filename,'_');
    expDate{i} = strcat(filenameSplit{2:4});
    expDateSimp{i} = strcat(filenameSplit{3},filenameSplit{4},filenameSplit{2}(3:4));
end

tempSessionNum = 0;
for i = 1:length(filenames)
    tic;
    
    filename = filenames(i).name;
    tempDate = expDate{i}; sameDateFlag = find(strcmp(tempDate,expDate));

    fileData = readtable([csvPath sep filename]);

    entryIdx = strcmp(fileData.Subject,'Entry');
    entryName = fileData.Var5(entryIdx); entryTime = fileData.Var2(entryIdx);

    inputIdx = strcmp(fileData.Subject,'Input');
    inputName = fileData.Var5(inputIdx); inputTime = fileData.Var2(inputIdx);
    var4Name = fileData.Var4(inputIdx); var4Time = fileData.Var2(inputIdx);


    if length(sameDateFlag)>1 % multiple files for one day
        tempEntryIdx = [tempEntryIdx ; entryIdx]; 
        tempEntryName = [tempEntryName ; entryName]; 
        tempEntryTime = [tempEntryTime ; entryTime + 100000*tempSessionNum ];

        tempInputIdx = [tempInputIdx ; inputIdx]; 
        tempInputName = [tempInputName ; inputName]; 
        tempInputTime = [tempInputTime ; inputTime + 100000*tempSessionNum ];
        tempSessionNum = tempSessionNum+1;

        tempVar4Name = [tempVar4Name ; var4Name]; 
        tempVar4Time = [tempVar4Time ; var4Time]; 
    else
        tempEntryIdx = entryIdx; tempEntryName = entryName; tempEntryTime = entryTime;
        tempInputIdx = inputIdx; tempInputName = inputName; tempInputTime = inputTime;
        tempVar4Name = var4Name; tempVar4Time = var4Time;
    end

    if length(sameDateFlag)==1 || (length(sameDateFlag)>1 && sameDateFlag(end) == i)
        trainDay = trainDay + 1;

        % STIMULUS category, 1 or 2
        trialData.stimulus{trainDay} = fn_removeNan(fn_multistrcmpCategory(tempEntryName, stateKey.stimulusKey))';
        stimFlag = logical(sum(~isnan(fn_multistrcmpCategory(tempEntryName, stateKey.stimulusKey)),1));
        % CONTEXT category, 1 - reinf, 2 - correction, 3 - probe
        trialData.context{trainDay} = fn_removeNan(fn_multistrcmpCategory(tempEntryName, stateKey.contextKey))';
        % RESPONSE type category, 0 - miss, 1 - correct, 2 - incorrect
        trialData.responseType{trainDay} = fn_removeNan(fn_multistrcmpCategory(tempEntryName, stateKey.responseKey))'-1;
        respFlag = logical(sum(~isnan(fn_multistrcmpCategory(tempEntryName, stateKey.responseKey)),1));
        % ACTION category, 0 - no action, 1 - left, 2 - right
        trialData.action{trainDay} = fn_removeNan(fn_multistrcmpCategory(tempEntryName, stateKey.actionKey))'-1;
        % REWARD category, 0 - no reward, 1 - reward
        trialData.reward{trainDay} = fn_removeNan(fn_multistrcmpCategory(tempEntryName, stateKey.rewardKey))'-1;

        % CORRECT FOR DAYS WHERE PROTOCOl TERMINATES BEFORE RESPONSE 
        if sum(stimFlag) > sum(respFlag)
            if (sum(stimFlag) - sum(respFlag)) == 1
                disp('Session ended prematurely')
                trialData.stimulus{trainDay}(end) = [];
                trialData.context{trainDay}(end) = [];
                tempIdx = find(stimFlag==1); 
                stimFlag(tempIdx(end)) = 0;
            end
        end
        
        if sum(stimFlag) < sum(respFlag)
            disp('Probe bugged out');
            tempRespIdx = find(respFlag==1); 
            tempRespIdx(1:sum(stimFlag)) = [];
            respFlag(tempRespIdx) = 0;
            
            trialData.responseType{trainDay}(sum(stimFlag)+1:end) = [];
            trialData.action{trainDay}(sum(stimFlag)+1:end) = [];
            trialData.reward{trainDay}(sum(stimFlag)+1:end) = [];
        end
        
        % STIMULUS time 
        trialData.stimulusTime{trainDay} = tempEntryTime(stimFlag);
        % RESPONSE time
        trialData.responseTime{trainDay} = tempEntryTime(respFlag);
        
        % LICK time
        if ~isempty(trialData.stimulusTime{trainDay})
            if strcmp(mouse,'zz113') || strcmp(mouse,'zz115')
                tempLickTime = tempInputTime(tempVar4Name == 37);
            else
                tempLickTime = tempInputTime(strcmp(tempInputName,'LickOn'));
            end
            tempStimTime = trialData.stimulusTime{trainDay};
            tempSessionIdxStim = find(diff(trialData.stimulusTime{trainDay})<0); 
            tempSessionIdxLick = find(diff(tempLickTime)<0);
            tempSessionIdxStim = [tempSessionIdxStim; length(tempStimTime)];
            tempSessionIdxLick = [tempSessionIdxLick; length(tempLickTime)];
            %{
            for j = 1:length(tempSessionIdxStim)-1
                tempIdx = tempSessionIdxStim(j)+1:tempSessionIdxStim(j+1);
                tempStimTime(tempIdx) = tempStimTime(tempIdx) + 100000*j;
                tempIdx = tempSessionIdxLick(j)+1:tempSessionIdxLick(j+1);
                tempLickTime(tempIdx) = tempLickTime(tempIdx) + 100000*j;
            end
            %}

            tempStimTime = [tempStimTime; tempStimTime(end)+100];

            [nLick,~,lickTrialIdx] = histcounts(tempLickTime,tempStimTime);
            disp(['This day have n=' int2str(length(sameDateFlag)) ' files, maxLick=' int2str(max(nLick))])
            trialData.nlick{trainDay} = nLick';
            lickTime = {}; lickTimeStr = {};
            for j = 1:length(trialData.stimulusTime{trainDay})
                lickTime{j} = tempLickTime(lickTrialIdx==j) - tempStimTime(j);
                lickTimeStr{j} = string(strjoin(arrayfun(@(x)(num2str(x,'%.3f')), lickTime{j}, 'UniformOutput', false)));
            end 
            try
                trialData.lickTime{trainDay} = fn_cell2matFillNan(lickTime)';
                if size(trialData.lickTime{trainDay},1) == 1
                    trialData.lickTime{trainDay} = trialData.lickTime{trainDay}'; end
            catch
                disp('uh oh')
            end
            trialData.lickTimeStr{trainDay} = lickTimeStr';
        end

        % DATE
        trialData.date{trainDay} = expDate{i};
        % TRAINING TYPE
        trialData.trainingType{trainDay} = dirName;
        
        % Wheel information
        filePath = [wheelTxtPath sep mouse '_' expDateSimp{i} '.txt'];
        wheelSoundOnCheckFlag = 0;
        disp('------------------------------------------------------------')
        if any(contains({'zz054','zz062','zz063','zz066','zz067','zz068','zz069'},mouse))
            if isfile(filePath)
                disp([mouse '_' expDateSimp{i} '; nTrial = ' int2str(length(trialData.action{trainDay})) '; ' dirName])
                if strcmp(dirName,'wheelTraining') % For wheel training Txt
                    wheelBout = readWheelTxt(filePath);
                    trialData.wheelSoundOn{trainDay} = {};
                    trialData.wheelPreSound{trainDay} = [];
                    trialData.wheelSoundOnCheckFlag{trainDay} = wheelSoundOnCheckFlag;
                else % For 2AFC Txt
                    [wheelSoundOn,wheelPreSoundTemp,hitCool] = read2AFCTxt(filePath);  
                    wheelPreSound = cellfun(@cell2mat,wheelPreSoundTemp,'UniformOutput',false);
                    wheelPreSound = cellfun(@length,wheelPreSound)';
                    hitCoul = trialData.action{trainDay} ~= 0;
                    if length(hitCool) == length(hitCoul)
                        sameFlag = hitCool==hitCoul'; 
                        disp(['Coulbourn & Coolterm SAME Trialnum (' int2str(length(hitCoul)) '); Hit Check = ' int2str(sum(sameFlag))]);
                        if sum(sameFlag)==length(hitCool); wheelSoundOnCheckFlag = 1; 
                        %elseif length(hitCool)- sum(sameFlag) < 5 && all(hitCoul(sameFlag==0)==0)
                        %    disp([int2str(length(hitCool)- sum(sameFlag)) ' trials have actions recorded as miss']);                        
                        %    wheelSoundOnCheckFlag = 1; 
                        elseif length(hitCool)- sum(sameFlag) < 10
                            wheelSoundOn(sameFlag==0) = cell(1,sum(sameFlag==0)); wheelSoundOnCheckFlag = 1;
                            wheelPreSound(sameFlag==0) = nan;
                            disp([int2str(length(hitCool)- sum(sameFlag)) ' trials wheel data discarded']); 
                        end


                    else
                        disp(['Coulbourn & Coolterm DIFFERENT Trialnum (' int2str(length(hitCoul)) ' ' int2str(length(hitCool)) ')'])
                        if length(hitCool) > length(hitCoul)
                            tempTrialDiff = length(hitCool) - length(hitCoul);
                            if all(cellfun(@isempty,wheelSoundOn(1:tempTrialDiff))) && ...
                                    all(cellfun(@isempty,wheelPreSoundTemp(1:tempTrialDiff)))
                                wheelSoundOn = wheelSoundOn(tempTrialDiff+1:end);
                                wheelPreSound = wheelPreSound(tempTrialDiff+1:end);
                                hitCool = hitCool(tempTrialDiff+1:end);
                                sameFlag = hitCool==hitCoul';
                                disp(['FIXED -- Trialnum different likely due to sound test before 2AFC. Hit Check = ' int2str(sum(sameFlag))])
                                if sum(sameFlag)==length(hitCool); wheelSoundOnCheckFlag = 1;
                                %elseif length(hitCool)- sum(sameFlag) < 5 && all(hitCoul(sameFlag==0)==0)
                                %    disp([int2str(length(hitCool)- sum(sameFlag)) ' trials have actions recorded as miss']);                        
                                %    wheelSoundOnCheckFlag = 1; 
                                elseif length(hitCool)- sum(sameFlag) < 10
                                    wheelSoundOn(sameFlag==0) = cell(1,sum(sameFlag==0)); wheelSoundOnCheckFlag = 1;
                                    wheelPreSound(sameFlag==0) = nan;
                                    disp([int2str(length(hitCool)- sum(sameFlag)) ' trials wheel data discarded']); 
                                end
                            end
                        end



                    end
                    if ~isempty(wheelSoundOn)
                        trialData.wheelSoundOn{trainDay} = fn_cell2matFillNan(wheelSoundOn);
                        trialData.wheelPreSound{trainDay} = wheelPreSound;
                    else 
                        trialData.wheelSoundOn{trainDay} ={}; trialData.wheelPreSound{trainDay} =[];
                    end
                    trialData.wheelSoundOnCheckFlag{trainDay} = wheelSoundOnCheckFlag;
                end
            else
                disp([mouse '_' expDateSimp{i} '; ' dirName '; Not wheel training or 2AFC session OR txt missing'])
                trialData.wheelSoundOn{trainDay} = {};trialData.wheelPreSound{trainDay} = [];
                trialData.wheelSoundOnCheckFlag{trainDay} = wheelSoundOnCheckFlag;
            end
        else
            if isfile(filePath)
                disp([mouse '_' expDateSimp{i} '; nTrial = ' int2str(length(trialData.action{trainDay})) '; ' dirName])
                %try
                    [trialData.wheelPos_aligned{trainDay},~] = fn_read2AFCTxtTimeStamp(filePath);
                %catch
                %    trialData.wheelPos_aligned{trainDay} = [];
                %    disp('TXT reading code failed')
                %end
            else
                disp([mouse '_' expDateSimp{i} '; ' dirName '; Not wheel training or 2AFC session OR txt missing'])
                trialData.wheelPos_aligned{trainDay} =[];
            end
        end
        % REACTION TIME
        try
            trialData.reactionTime{trainDay} = trialData.responseTime{trainDay}-trialData.stimulusTime{trainDay};
            disp(['Mean reaction time = ' num2str(mean(trialData.reactionTime{trainDay}))])
        catch
           error('Session ended prematurely'); 
        end
        % RESET variables
        tempEntryIdx = []; tempEntryName = []; tempEntryTime = []; 
        tempInputIdx = []; tempInputName = []; tempInputTime = []; 
        tempSessionNum = 0;
        tempVar4Name = []; tempVar4Time = [];
    end
    toc;
end


end


function wheelBout = readWheelTxt(txtFilename)
global sep mouse;

wheelBout = {}; 
if isfile(txtFilename)
    txtData = splitlines(fileread(txtFilename));
    for i = 1:length(txtData)
        if ~isempty(txtData{i}) && ~isletter(txtData{i}(1)) 
            % Wheel position if the first char is not letter 
            txtDataSplit = strsplit(txtData{i});
            txtDataSplitMat = cell2mat(cellfun(@str2double,txtDataSplit,'UniformOutput',false));
            txtDataSplitMat(isnan(txtDataSplitMat)) = [];
            wheelBout{end+1} = wheelBout;
        end

    end

else
    error(['Error - txtFilename ' txtFilename ' does not exist!']);
end

end


function [wheelSoundOn,wheelPreSound,hit] = read2AFCTxt(txtFilename)
global sep mouse;

nTrial = 0; soundOnFlag = false; hit = [];
wheelSoundOn = {}; wheelPreSound = {};  
tempSound = {};  tempPreSound = {}; 
if isfile(txtFilename)
    txtData = splitlines(fileread(txtFilename));
    for i = 1:length(txtData)                  
        if ~isempty(txtData{i}) && ~isletter(txtData{i}(1)) % Wheel position if the first char is not letter
            txtDataSplit = strsplit(txtData{i});
            txtDataSplitMat = cell2mat(cellfun(@str2double,txtDataSplit,'UniformOutput',false));
            if sum(isnan(txtDataSplitMat)) ~=0 
                % Get rid of all entries after the first nan (i.e. letter)
                % Due to the exception of 'Sound 1 on!'
                firstNan = fn_findFirst(isnan(txtDataSplitMat));
                txtDataSplitMat(firstNan:end) = [];
            end
            if soundOnFlag
                tempSound{end+1} = txtDataSplitMat;
            else
                tempPreSound{end+1} = txtDataSplitMat;
            end      
        end
        
        
        % Since 'Sonud on' and 'Sound off' always at the end of a line,
        % record wheel first, then check if song is on or off
        if contains(txtData{i},'Sound')
            if contains(txtData{i},'on')
                soundOnFlag = true; nTrial = nTrial+1;
                wheelPreSound{nTrial} = tempPreSound; 
                tempPreSound = {}; tempSound = {};
            elseif contains(txtData{i},'off')
                if nTrial == 0; nTrial = nTrial + 1; end % In case of sound-on not recorde in trial 1
                soundOnFlag = false; hit(nTrial) = 0;
                if ~isempty(tempSound); wheelSoundOn{nTrial} = tempSound{1}; 
                else; wheelSoundOn{nTrial} = []; end
                if length(tempSound) >1
                    disp('WARNING -- multiple movement in sound on'); 
                end
                tempSound = {}; tempPreSound = {};
            end
        elseif contains(txtData{i},'choice reached')
            % Choice reached always follow sound off; only record action
            hit(nTrial) = 1;
        end

    end

    
else
    error(['Error - txtFilename ' txtFilename ' does not exist!']);
end


end

