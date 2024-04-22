function [transitionData, transitionID] = fn_getTransition(block)
    transitionStart = find(diff(block.stateFlag)~=0);
    blockCount = [0 0 0]; 
    transitionData = zeros(length(transitionStart),4);
    transitionID = cell(length(transitionStart),1);
    for i = 1:length(transitionStart)
        prevBlock = block.stateFlag(transitionStart(i)); 
        nextBlock = block.stateFlag(transitionStart(i)+1);   
        [trialMat, transID] = updateTransitionData(blockCount,block,prevBlock,nextBlock);
        blockCount = updateBlockCount(blockCount,prevBlock);
        transitionData(i,:) = trialMat;
        transitionID{i} = transID;
    end
    
    
    %temp = {idxA.start,idxA.end,idxB.start,idxB.end};
    %[A2B] = getIdx(idxA.start,idxA.end,idxB.start,idxB.end); % A to B
    %[B2A] = getIdx(idxB.start,idxB.end,idxA.start,idxA.end); % B to A

    function [A2B] = getIdx(startA,endA,startB,endB)
        tempDiff = repmat(startB',length(endA),1)-repmat(endA,1,length(startB)); 
        A2B = nan(length(endA),4); 
        for i = 1:length(endA)
            tempIdx = find(tempDiff(i,:) > 0); 
            if ~isempty(tempIdx)
                tempIdx = tempIdx(1);
                A2Bdiff(i) = tempDiff(i,tempIdx);
                A2B(i,1) = startA(i); A2B(i,2) = endA(i); 
                A2B(i,3) = startB(tempIdx); A2B(i,4) = endB(tempIdx);
            end
        end
    end
end

function [trialMat, transID] = updateTransitionData(blockCount,blocks,prevBlock,nextBlock)
    trialMat = zeros(1,4);
    
    if prevBlock==1
        trialMat(1) =  blocks.blockL.start(blockCount(1)+1);
        trialMat(2) =  blocks.blockL.end(blockCount(1)+1);
        letter1 = 'L';
    elseif prevBlock==-1
        trialMat(1) =  blocks.blockR.start(blockCount(2)+1);
        trialMat(2) =  blocks.blockR.end(blockCount(2)+1);
        letter1 = 'R';
    elseif prevBlock==0
        trialMat(1) =  blocks.blockU.start(blockCount(3)+1);
        trialMat(2) =  blocks.blockU.end(blockCount(3)+1);
        letter1 = 'U';
    end

    if nextBlock==1
        trialMat(3) =  blocks.blockL.start(blockCount(1)+1);
        trialMat(4) =  blocks.blockL.end(blockCount(1)+1);
        letter2 = 'L';
    elseif nextBlock==-1
        trialMat(3) =  blocks.blockR.start(blockCount(2)+1);
        trialMat(4) =  blocks.blockR.end(blockCount(2)+1);
        letter2 = 'R';
    elseif nextBlock==0
        trialMat(3) =  blocks.blockU.start(blockCount(3)+1);
        trialMat(4) =  blocks.blockU.end(blockCount(3)+1);
        letter2 = 'U';
    end
    transID = strjoin({letter1, letter2},'2');


end

function blockCount = updateBlockCount(blockCount,prevBlock)

    if prevBlock==1; blockCount(1) =  blockCount(1) + 1; 
    elseif prevBlock==-1; blockCount(2) =  blockCount(2) + 1;
    elseif prevBlock==0; blockCount(3) =  blockCount(3) + 1;
    end

end

    