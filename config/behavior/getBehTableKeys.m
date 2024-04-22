function tableKeys = getBehTableKeys(behavTpe)
switch behavTpe
    case 'wheel2AFC'
    tableKeys = {'day','stimulus','context','responseType','action','reward','stimulusTime',...
        'responseTime','reactionTime','nlick','lickTimeStr','date','trainingType'};
end
    
end