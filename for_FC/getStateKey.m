function stateKey = getStateKey()
stateKey = struct();

stateKey.stimulusKey = {{'CueLeft'},{'CueRight'}};

stateKey.contextKey = {{'CueLeft'},{'CueRight'}};

stateKey.actionKey = {{'MissLeft','MissRight'},{'HitLeft'},{'HitRight'}};

end