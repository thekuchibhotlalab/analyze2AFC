function  readCoulbourn(mousename)
mouse = mousename;
dataPath = 'C:\Users\zzhu34\Documents\gitRep\analyze2AFCObj\for_FC';
writeDataPath = [dataPath filesep 'trialData' filesep]; mkdir(writeDataPath);

readDataPath = 'C:\Users\zzhu34\Documents\gitRep\analyze2AFCObj\for_FC';
mousePath = [readDataPath filesep mouse filesep];

%% read animal name and date for each file
dirName = fn_readDir(mousePath);
trialData = {};
for i = 1:length(dirName)
   trialData{i} = readFiles([mousePath filesep dirName{i}],dirName{i}); 
end
if length(trialData) > 1; trialData = fn_catStructField(2, trialData{:}); 
else; trialData = trialData{1};end
trialData = fn_sortStructByFieldKey(trialData,'date');

save([writeDataPath filesep mouse '.mat'],'trialData');
end
%% All the functions

function trialData = readFiles(csvPath,dirName)
 
filenames = dir([csvPath filesep '*.csv']);
stateKey = getStateKey();
trialName = []; responseName = [];

trialData = struct();
trainDay = 0; tempDate = [];
tempEntryIdx = []; tempEntryName = []; tempEntryTime = []; 
expDate = {}; expDateSimp = {};
for i = 1:length(filenames)
    filename = filenames(i).name;
    filenameSplit = strsplit(filename,'_');
    expDate{i} = strcat(filenameSplit{2:4});
    expDateSimp{i} = strcat(filenameSplit{3},filenameSplit{4},filenameSplit{2}(3:4));
end

for i = 1:length(filenames)
    tic;
    
    filename = filenames(i).name;
    tempDate = expDate{i}; sameDateFlag = find(strcmp(tempDate,expDate));

    fileData = readtable([csvPath filesep filename]);

    entryIdx = strcmp(fileData.Subject,'Entry');
    entryName = fileData.Var5(entryIdx); entryTime = fileData.Var2(entryIdx);



    if length(sameDateFlag)>1 % multiple files for one day
        tempEntryIdx = [tempEntryIdx ; entryIdx]; 
        tempEntryName = [tempEntryName ; entryName]; 
        tempEntryTime = [tempEntryTime ; entryTime];
    else
        tempEntryIdx = entryIdx; tempEntryName = entryName; tempEntryTime = entryTime;
    end

    if length(sameDateFlag)==1 || (length(sameDateFlag)>1 && sameDateFlag(end) == i)
        trainDay = trainDay + 1;

        % STIMULUS category, 1 or 2
        trialData.stimulus{trainDay} = fn_removeNan(fn_multistrcmpCategory(tempEntryName, stateKey.stimulusKey))';
        stimFlag = logical(sum(~isnan(fn_multistrcmpCategory(tempEntryName, stateKey.stimulusKey)),1));
        % CONTEXT category, 1 - reinf, 2 - correction, 3 - probe
        trialData.context{trainDay} = fn_removeNan(fn_multistrcmpCategory(tempEntryName, stateKey.contextKey))';
        
        % ACTION category, 0 - no action, 1 - left, 2 - right
        trialData.action{trainDay} = fn_removeNan(fn_multistrcmpCategory(tempEntryName, stateKey.actionKey))'-1;
        actionFlag = logical(sum(~isnan(fn_multistrcmpCategory(tempEntryName, stateKey.actionKey)),1));
        % STIMULUS time 
        trialData.stimulusTime{trainDay} = tempEntryTime(stimFlag);
        % RESPONSE time
        trialData.responseTime{trainDay} = tempEntryTime(actionFlag);
        
        % DATE
        trialData.date{trainDay} = expDate{i};
        % TRAINING TYPE
        trialData.trainingType{trainDay} = dirName;
        
    end
    toc;
end
end