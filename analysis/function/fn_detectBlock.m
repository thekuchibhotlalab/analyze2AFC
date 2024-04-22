function [stateFlag, highBlock, lowBlock, normBlock] = fn_detectBlock(mat,baseline, threshold,varargin)

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('blockIntervalThre', 5);
p.addParameter('blockLenThre', 5);
p.parse(varargin{:});

stateFlag = zeros(size(mat));
stateFlag (mat >= (baseline+threshold)) = 1; stateFlag (mat <= (baseline-threshold)) = -1;
% first correct out small blocks of high and low by blockIntervalThre
stateFlag = correctSmallBlocks(stateFlag,1,p.Results.blockLenThre);
stateFlag = correctSmallBlocks(stateFlag,-1,p.Results.blockLenThre);
% then detect high and low blocks, and combine blocks by blockLenThre
stateFlag = combineSmallBlocks(stateFlag,1,p.Results.blockIntervalThre);
stateFlag = combineSmallBlocks(stateFlag,-1,p.Results.blockIntervalThre);
% finally output (corrected) stateFlag and blocks

[highBlock.start, highBlock.end, highBlock.len] = fn_findBlock(stateFlag,1,1);
[lowBlock.start, lowBlock.end, lowBlock.len] = fn_findBlock(stateFlag,-1,1);
[normBlock.start, normBlock.end, normBlock.len] = fn_findBlock(stateFlag,0,1);

end


function stateFlag = correctSmallBlocks(stateFlag,element,blockLenThre)
    [startIdx,endIdx,len] = fn_findBlock(stateFlag,element,1);
    blockCorrIdx = find(len<=blockLenThre);
    for i = 1:length(blockCorrIdx); stateFlag(startIdx(blockCorrIdx(i)):endIdx(blockCorrIdx(i))) = 0; end
end

function stateFlag = combineSmallBlocks(stateFlag,element,blockIntervalThre)
    [startIdx,endIdx,~] = fn_findBlock(stateFlag,element,blockIntervalThre);
    for i = 1:length(startIdx); stateFlag(startIdx(i):endIdx(i)) = element; end
end
