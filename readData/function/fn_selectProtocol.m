function trialData = fn_selectProtocol(dataPath, mouse, selectProtocol)
global sep;
if isempty(sep); sep = '/'; end
if ~iscell(selectProtocol); selectProtocol = {selectProtocol}; end
load([dataPath sep mouse '.mat'],'trialData');
trialData = fn_readStructByFieldKey(trialData,'trainingType',selectProtocol);

end