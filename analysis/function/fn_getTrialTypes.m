function idxData = fn_getTrialTypes(stim, resp)

    
corrL = sum((stim==1) & (resp==1));
incorrL = sum((stim==1) & (resp==2));
missL = sum((stim==1) & (resp==0));

corrR = sum((stim==2) & (resp==1));
incorrR = sum((stim==2) & (resp==2));
missR = sum((stim==2) & (resp==0));

idxData = [corrL,incorrL,missL,corrR,incorrR,missR];


end