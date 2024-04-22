function allMouse = fn_getObjPT_bin30(mice,opsParam)
if nargin == 1; opsParam = struct(); opsParam.biasBlockType = 'AUC'; end
allMouse = {};
for i = 1:length(mice)
    
    ops = opsSharedArg(opsParam);
    ops = opsIndivArg(ops,mice{i});
    
    load([ops.datapath mice{i} '.mat']);

    tempObj = wheel2AFC(trialData,'ops',ops,'mouse',mice{i} );

    tempObj = tempObj.removeMiss(); 
    tempObj = tempObj.getAcc();
    tempObj = tempObj.calculateProbe();
    tempObj = tempObj.getPsyTrack();
    if strcmp(ops.biasBlockType,'threshold')
        tempObj = tempObj.getBiasBlock(0.40);
    elseif strcmp(ops.biasBlockType,'AUC')
        tempObj = tempObj.getBiasBlock2();
    end
    allMouse{i} = tempObj;
end
    
    
end


function ops = opsSharedArg(ops)
    
% ------------------ required arguments ------------------
ops.datapath = 'C:\Users\zzhu34\Documents\tempdata\octoData\trialData\';
ops.behavType = 'wheel2AFC';
ops.selectProtocol = {'puretone'};
ops.preprocessedInput = false;
ops.psytrackPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData\psyTrackFit\';
ops.correctWheelTrial = true;
ops.learningCurveBin = 30; 

end

function ops = opsIndivArg(ops,mouse)
    
switch mouse
    case 'zz054'; ops.startDay = []; ops.endDay = 20210610;
    case 'zz062'; ops.startDay = []; ops.endDay = 20210621;
    case 'zz063'; ops.startDay = []; ops.endDay = 20210611;
    case 'zz066'; ops.startDay = []; ops.endDay = 20210619;
    case 'zz067'; ops.startDay = []; ops.endDay = 20210706;
    case 'zz068'; ops.startDay = []; ops.endDay = 20210622;
    case 'zz069'; ops.startDay = []; ops.endDay = 20210629;
    case 'zz107'; ops.startDay = []; ops.endDay = 20220127;
    case 'zz109'; ops.startDay = []; ops.endDay = 20220127;
    case 'zz111'; ops.startDay = []; ops.endDay = 20220312;
    case 'zz112'; ops.startDay = []; ops.endDay = 20220312;
    case 'zz113'; ops.startDay = []; ops.endDay = 20220309;
    case 'zz115'; ops.startDay = []; ops.endDay = 20220312;
end


    
end