function S = fn_preprocessInput(S,ops)
%myFun - Description
%
% Syntax: output = myFun(input)
%
% Long description


if isfield(ops,'selectProtocol') && ~isempty(ops.selectProtocol); S = selectProtocol(S,ops.selectProtocol); end
if isfield(ops,'startDay') && ~isempty(ops.startDay); S = selectStartDay(S,ops.startDay); end
if isfield(ops,'endDay') && ~isempty(ops.endDay); S = selectEndDay(S,ops.endDay); end

for i = 1:length(S.stimulus)
    S.day{i} = ones(size(S.stimulus{i})) * i;
end

tempField  = fieldnames(S);
for i = 1:length(tempField)
    if ischar(S.(tempField{i}){1}) || isstring(S.(tempField{i}){1})
        tempSize = cellfun(@(x)size(x,1),S.stimulus,'UniformOutput',false);
        S.(tempField{i}) = cellfun(@(x,y)repmat(string(x),y,1),(S.(tempField{i})),tempSize,'UniformOutput',false);
    end
end
S = concatenateDay(S,ops);

end

function S = selectProtocol(S,protocol)
    S = selectIdx(S,fn_multistrcmp(S.trainingType,protocol));
end

function S = selectStartDay(S,startDay)
    tempDate = cellfun(@str2double,S.date); S = selectIdx(S,tempDate>=startDay);
end

function S = selectEndDay(S,endDay)
    tempDate = cellfun(@str2double,S.date); S = selectIdx(S,tempDate<=endDay);
end

function S = concatenateDay(S,ops)
    S = fn_structfun(@(x)fn_cell2mat(x,1),S,'funcName', 'concatenateDay','verbose',ops.verbose);
end

function S = selectIdx(S,dayFlag)
    S = fn_readStructBySelection(S,'selectFlag',dayFlag);
end