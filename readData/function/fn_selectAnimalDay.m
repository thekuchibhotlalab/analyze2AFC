
function trialData = fn_selectAnimalDay(trialData,mouse)

switch mouse
    case 'zz054'
        maxDay = 20210610;
    case 'zz062'
        maxDay = 20210621;
    case 'zz063'
        maxDay = 20210611;
    
    case 'zz066'
        maxDay = 20210619;
    case 'zz067'
        maxDay = 20210706;
    case 'zz068'
        maxDay = 20210622;
    case 'zz069'
        maxDay = 20210629;
    case 'zz097'
        maxDay = 20220125;
    case 'zz098'
        maxDay = 20220201;
    case 'zz099'
        maxDay = 20220121;
    case 'zz100'
        maxDay = 20220121;
    case 'zz101'
        maxDay = 20220121;
    case 'zz102'
        maxDay = 20220121;
    case 'zz105'
        maxDay = 20220121;
    case 'zz107'
        maxDay = 20220127;
    case 'zz109'
        maxDay = 20220201;
    case 'zz111'
        maxDay = 20220312;
    case 'zz112'
        maxDay = 20220312;
    case 'zz113'
        maxDay = 20220312;
    case 'zz115'
        maxDay = 20220312;
            
end

tempDate = cellfun(@str2double,trialData.date);
trialData = fn_readStructByFlag(trialData,tempDate<=maxDay);

end