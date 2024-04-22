function stimProb = fn_getStimProb(mouseObj)
    stimProb = smoothdata(mouseObj.behav.stimulus==1,'movmean',...
        mouseObj.ops.learningCurveBin);
end