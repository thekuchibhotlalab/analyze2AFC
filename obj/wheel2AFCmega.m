classdef wheel2AFCmega 
    %WHEEL2AFCMEGA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mouseCell
        nMouse 
        nTrials
        mouseName
        
    end
    
    methods
        function obj = wheel2AFCmega(objCell)
            %WHEEL2AFCMEGA Construct an instance of this class
            
            obj.mouseCell = objCell; obj.nMouse = length(objCell);
            obj.nTrials = cellfun(@(x)(size(x.behav,1)),objCell);
            obj.mouseName = cellfun(@(x)(x.ops.mouse),objCell,'UniformOutput',false);
            %   Detailed explanation goes here
            %{
            objPropName = properties('wheel2AFCmega');
            for i = 1:length(objPropName)
                if any(strcmp(objPropName{i}, properties(objCell{1})))
                    obj.(objPropName{i}) = combineObj(objCell,objPropName{i});
                end
            end
            obj.nMouse = length(obj.(objPropName{1}));
            function C =  combineObj(objCell,prop)
                C = cell(1,length(objCell)); for j = 1:length(objCell); C{j} = objCell{j}.(prop); end
            end
            %}
        end
        
        function outCell = getProp(obj,prop,varargin)
            p = inputParser; p.KeepUnmatched = true; arg = {'matFlag','field','idx';false,'',[]}; 
            cellfun(@(x,y)(p.addParameter(x,y)),arg(1,:),arg(2,:)); p.parse(varargin{:});
            %if ~contains(propName, obj.behav{1}.Properties.VariableNames); msgbox('Input property not in wheel2AFC object.','ERROR'); return; end
            if ~isempty(p.Results.field); outCell = cellfun(@(x)(x.(prop).(p.Results.field)),obj.mouseCell,'UniformOutput',false); 
            else; outCell = cellfun(@(x)(x.(prop)),obj.mouseCell,'UniformOutput',false); 
            end

            if ~isempty(p.Results.idx); outCell = outCell(p.Results.idx); end
            if p.Results.matFlag; outCell = fn_cell2matFillNan(outCell); end
        end
        
        function reinf = loadReinf(obj)
            reinf.acc = obj.getProp('behav','field','acc','matFlag',true);
            reinf.bias = abs(obj.getProp('behav','field','bias','matFlag',true));
            reinf.biasDir = obj.getProp('behav','field','bias','matFlag',true);
            reinf.dp = abs(obj.getProp('behav','field','dprime','matFlag',true));
            reinf.actionRate = abs(obj.getProp('behav','field','actionRate','matFlag',true));
            reinf.reactionTime = abs(obj.getProp('behav','field','reactionTime','matFlag',true));
            %acc.probe = obj.getProp('probe','field','probeAcc','matFlag',true);
            %bias.probe = obj.getProp('probe','field','probeBias','matFlag',true);
            %probeIdx = obj.getProp('probe','field','probeIdxBin','matFlag',true);
            %acc.bef = obj.getProp('probe','field','befAcc','matFlag',true);
            %bias.bef = obj.getProp('probe','field','befBias','matFlag',true);
            %acc.aft = obj.getProp('probe','field','aftAcc','matFlag',true);
            %bias.aft = obj.getProp('probe','field','aftBias','matFlag',true); 
        end

        function tblCell = loadAnimalBehav(obj,nameCell,varargin)
            p = inputParser; p.KeepUnmatched = true; arg = {'idx';[]}; 
            cellfun(@(x,y)(p.addParameter(x,y)),arg(1,:),arg(2,:)); p.parse(varargin{:});
            
            tblCell = cellfun(@(x)(fn_readTableVar(x.behav,nameCell)),obj.mouseCell,'UniformOutput',false);
            if ~isempty(p.Results.idx); tblCell = tblCell(p.Results.idx); end
        end
        
        function probeBin = loadProbeBin(obj)
            outCell = obj.objFun('binProbe',{});
            probeBin = fn_catStructField(2, outCell{1}{:});
        end
        
        function outCell = objFun(obj,funcName,funcInputs,varargin)
            p = inputParser; p.KeepUnmatched = true; arg = {'cellfun','fun';[],[]}; 
            cellfun(@(x,y)(p.addParameter(x,y)),arg(1,:),arg(2,:)); p.parse(varargin{:});
            
            
            fhandle = @(y)(cellfun(@(x)(x.(funcName)(funcInputs{:})),y,'UniformOutput',false)); 
            try nargs = nargout(['wheel2AFC>wheel2AFC.' funcName] );
                catch; nargs = nargout(['twoChoiceTask>twoChoiceTask.' funcName] ); end
            outCell = fn_getFnOut(fhandle,{obj.mouseCell},'nargs',nargs);
            
            if ~isempty(p.Results.cellfun)
                for i = 1:length(outCell)
                    outCell{i} = cellfun(p.Results.cellfun,outCell{i},'UniformOutput',false);
                end
            end
            
            if ~isempty(p.Results.fun)
                for i = 1:length(outCell)
                    outCell{i} = p.Results.fun(outCell{i});
                end
                
            end
            
        end
    end
end

